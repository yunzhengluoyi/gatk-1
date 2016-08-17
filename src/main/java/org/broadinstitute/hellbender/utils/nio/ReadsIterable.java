package org.broadinstitute.hellbender.utils.nio;

import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.seekablestream.ByteArraySeekableStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.nio.ChannelAsSeekableStream;
import org.broadinstitute.hellbender.utils.nio.SeekableByteChannelPrefetcher;

import java.io.IOException;
import java.io.Serializable;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Iterator;
import java.util.List;
import java.util.NoSuchElementException;

/**
 * ReadsIterable gives you all the reads for a given genomic interval.
 *
 * QueryInterval + header --> iterable SAMRecords
 */
public class ReadsIterable implements Iterable<SAMRecord>, Serializable {

    private static final long serialVersionUID = 1L;
    final String path;
    final byte[] index;
    final QueryInterval interval;
    final boolean removeHeader = true;

    class ReadsIterator implements Iterator<SAMRecord> {
        final int BUFSIZE = 200 * 1024 * 1024;
        private SamReader bam;
        private long count = 0;
        private SAMRecordIterator query;
        private SAMRecord nextRecord = null;
        private boolean done = false;

        public ReadsIterator() throws IOException {
            List<SAMRecord> ret = new ArrayList<>();
            Path fpath = BucketUtils.getPathOnGcs(path);
            byte[] data = index;
            SeekableStream indexInMemory = new ByteArraySeekableStream(data);
            SeekableByteChannelPrefetcher chan2 = new SeekableByteChannelPrefetcher(Files.newByteChannel(fpath), BUFSIZE);
            ChannelAsSeekableStream bamOverNIO = new ChannelAsSeekableStream(chan2, path.toString());
            bam = SamReaderFactory.makeDefault()
                    .validationStringency(ValidationStringency.LENIENT)
                    .enable(SamReaderFactory.Option.CACHE_FILE_BASED_INDEXES)
                    .open(SamInputResource.of(bamOverNIO).index(indexInMemory));

            QueryInterval[] array = new QueryInterval[1];
            array[0] = interval;
            query = bam.query(array, false);
        }

        /**
         * Returns {@code true} if the iteration has more elements.
         * (In other words, returns {@code true} if {@link #next} would
         * return an element rather than throwing an exception.)
         *
         * @return {@code true} if the iteration has more elements
         */
        @Override
        public boolean hasNext() {
            if (done) return false;

            if (nextRecord!=null) return true;

            nextRecord = fetchRecord();

            boolean ret = (nextRecord != null);
            if (!ret) {
                done = true;
                try {
                    query.close();
                    query = null;
                    bam.close();
                    bam = null;
                } catch (IOException e) {
                    throw new RuntimeException(e);
                }
            }
            return ret;
        }

        /**
         * Returns the next element in the iteration.
         *
         * @return the next element in the iteration
         * @throws NoSuchElementException if the iteration has no more elements
         */
        @Override
        public SAMRecord next() {
            if (!hasNext()) throw new NoSuchElementException();
            SAMRecord ret = nextRecord;
            nextRecord = null;
            count++;
            return ret;
        }

        private SAMRecord fetchRecord() {
            while (query.hasNext()) {
                SAMRecord sr = query.next();
                int start = sr.getAlignmentStart();
                if (start >= interval.start && start <= interval.end) {
                    // read starts in the interval
                    if (removeHeader) {
                        sr.setHeader(null);
                    }
                    count++;
                    return sr;
                }
            }
            return null;
        }
    }

    public ReadsIterable(String path, byte[] index, QueryInterval in) {
        this.path = path;
        this.index = index;
        this.interval = in;
    }

    public Iterator<SAMRecord> iterator() {
        try {
            return new ReadsIterator();
        } catch (IOException x) {
            throw new RuntimeException(x);
        }
    }

}
