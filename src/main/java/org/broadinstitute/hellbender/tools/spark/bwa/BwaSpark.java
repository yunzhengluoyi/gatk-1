package org.broadinstitute.hellbender.tools.spark.bwa;

import htsjdk.samtools.SAMFileHeader;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.cmdline.Argument;
import org.broadinstitute.hellbender.cmdline.ArgumentCollection;
import org.broadinstitute.hellbender.cmdline.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.SparkProgramGroup;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.engine.filters.ReadFilterLibrary;
import org.broadinstitute.hellbender.engine.spark.GATKSparkTool;
import org.broadinstitute.hellbender.engine.spark.datasources.ReadsSparkSink;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadsWriteFormat;

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

@CommandLineProgramProperties(summary = "Runs BWA",
        oneLineSummary = "BWA on Spark",
        programGroup = SparkProgramGroup.class)
public final class BwaSpark extends GATKSparkTool {

    private static final long serialVersionUID = 1L;

    @Argument(doc = "the output bam", shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME,
            fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, optional = false)
    private String output;

    @ArgumentCollection
    private final BwaArgumentCollection bwaArgs = new BwaArgumentCollection();

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public boolean requiresReads() {
        return true;
    }

    @Override
    protected void runTool(final JavaSparkContext ctx) {
        final JavaRDD<GATKRead> unalignedReads = getReads();
        final String referenceFileName = referenceArguments.getReferenceFileName();
        final BwaSparkEngine engine = new BwaSparkEngine(bwaArgs.numThreads, bwaArgs.fixedChunkSize, referenceFileName);
        final SAMFileHeader readsHeader = engine.makeHeaderForOutput(getHeaderForReads(), getReferenceSequenceDictionary());
        final JavaRDD<GATKRead> reads = engine.alignWithBWA(ctx, unalignedReads, readsHeader);

        try {
            ReadsSparkSink.writeReads(ctx, output, null, reads, readsHeader, shardedOutput ? ReadsWriteFormat.SHARDED : ReadsWriteFormat.SINGLE);
        } catch (final IOException e) {
            throw new GATKException("Unable to write bam",e);
        }
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        final List<ReadFilter> filters = new ArrayList<>();
        filters.add(ReadFilterLibrary.HAS_MATCHING_BASES_AND_QUALS);
        filters.add(ReadFilterLibrary.SEQ_IS_STORED);
        return filters;
    }
}
