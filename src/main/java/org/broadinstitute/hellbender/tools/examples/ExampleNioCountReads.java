package org.broadinstitute.hellbender.tools.examples;

import com.google.cloud.storage.contrib.nio.CloudStorageFileSystem;
import com.google.common.base.Stopwatch;
import htsjdk.samtools.BAMFileSpan;
import htsjdk.samtools.BAMIndex;
import htsjdk.samtools.Chunk;
import htsjdk.samtools.QueryInterval;
import htsjdk.samtools.SAMFileHeader;
import htsjdk.samtools.SAMRecord;
import htsjdk.samtools.SAMRecordIterator;
import htsjdk.samtools.SAMSequenceRecord;
import htsjdk.samtools.SamInputResource;
import htsjdk.samtools.SamReader;
import htsjdk.samtools.SamReaderFactory;
import htsjdk.samtools.ValidationStringency;
import htsjdk.samtools.seekablestream.ByteArraySeekableStream;
import htsjdk.samtools.seekablestream.SeekableStream;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.ReadProgramGroup;
import org.broadinstitute.hellbender.engine.spark.SparkCommandLineProgram;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.nio.ChannelAsSeekableStream;
import org.broadinstitute.hellbender.utils.nio.NioBam;
import org.broadinstitute.hellbender.utils.nio.ReadsIterable;
import org.broadinstitute.hellbender.utils.nio.SeekableByteChannelPrefetcher;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.Path;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Calendar;
import java.util.List;
import java.util.concurrent.TimeUnit;

/**
 * Example of how to use Spark on Google Cloud Storage directly, without using the GCS Hadoop Connector.
 */
@CommandLineProgramProperties(
    summary = "Example of how to use Spark on Google Cloud Storage directly, without using the GCS Hadoop Connector",
    oneLineSummary = "Example of how to use Spark on Google Cloud Storage directly, without using the GCS Hadoop Connector",
    programGroup = ReadProgramGroup.class
)
public class ExampleNioCountReads extends SparkCommandLineProgram {
    private static final long serialVersionUID = 1L;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "Output file (if not provided, defaults to STDOUT)", common = false, optional = true)
    private File OUTPUT_FILE = null;

    @Argument(fullName = "inputPath", shortName = "P", doc = "Input path (eg. gs://foo/bar.bam)", optional = false)
    private String path = null;

    // Typically set to number of executors times number of cores per executor.
    @Argument(fullName = "parts", doc = "number of partitions", optional = false)
    private int parts = 3;

    private PrintStream outputStream;

    private void countReads(JavaSparkContext ctx) throws IOException {

        try {
            outputStream = OUTPUT_FILE != null ? new PrintStream(OUTPUT_FILE) : System.out;
        }
        catch ( FileNotFoundException e ) {
            throw new UserException.CouldNotReadInputFile(OUTPUT_FILE, e);
        }

        NioBam input = new NioBam(path, path + ".bai");
        long readCount = input.getReads(ctx, parts).count();
        outputStream.println("Number of reads: " + readCount);
    }

    /**
     * Runs the pipeline.
     *
     * @param ctx
     */
    @Override
    protected void runPipeline(JavaSparkContext ctx) {
        try {
            countReads(ctx);
        } catch (IOException x) {
            System.err.println(x);
        }
    }
}
