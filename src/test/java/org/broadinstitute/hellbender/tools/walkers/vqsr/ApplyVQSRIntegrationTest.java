package org.broadinstitute.hellbender.tools.walkers.vqsr;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.engine.FeatureDataSource;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.variant.GATKVariant;
import org.broadinstitute.hellbender.utils.variant.VariantContextVariantAdapter;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

public class ApplyVQSRIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return ApplyVQSR.class.getSimpleName();
    }

    public String getToolTestDataDir(){
        return publicTestDir + "org/broadinstitute/hellbender/tools/VQSR/";
    }

    public String getLargeVQSRTestDataDir(){
        return largeFileTestDir + "VQSR/";
    }

    @Test
    public void testApplySNPRecalibration() throws IOException {
        final String inputFile = getLargeVQSRTestDataDir() + "phase1.projectConsensus.chr20.1M-10M.raw.snps.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " -L 20:1,000,000-10,000,000" +
                    " --variant " + inputFile +
                    " --lenient" +
                    " --output %s" +
                    " -mode SNP" +
                    " -tranchesFile " + getLargeVQSRTestDataDir() + "expected/SNPTranches.txt" +
                    // pass in the tranche file to match GATK3; though without a TS_FILTER_LEVEL
                    // arg they aren't used
                    " -recalFile " + getLargeVQSRTestDataDir() + "snpRecal.vcf",
                Arrays.asList(getLargeVQSRTestDataDir() + "expected/snpApplyResult.vcf"));
        spec.executeTest("testApplyRecalibrationSNP", this);
    }

    @Test
    public void testApplyIndelRecalibrationIndel() throws IOException {
        final String inputFile = getLargeVQSRTestDataDir() + "combined.phase1.chr20.raw.indels.filtered.sites.1M-10M.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " -L 20:1,000,000-10,000,000" +
                    " -mode INDEL" +
                    " --lenient" +
                    " --variant " + inputFile +
                    " --output %s" +
                    // pass in the tranche file to match GATK3; though without a TS_FILTER_LEVEL
                    // arg they aren't used
                    " -tranchesFile " + getLargeVQSRTestDataDir() + "expected/indelTranches.txt" +
                    " -recalFile " + getLargeVQSRTestDataDir() + "indelRecal.vcf",
                Arrays.asList(getLargeVQSRTestDataDir() + "expected/indelApplyResult.vcf"));
        spec.executeTest("testApplyRecalibrationIndel", this);
    }

    @Test
    public void testApplyRecalibrationSnpAndIndelTogether() throws IOException {
        final IntegrationTestSpec spec = new IntegrationTestSpec(
                    " -L 20:1000100-1000500" +
                    " -mode BOTH" +
                    " --variant " + getToolTestDataDir() + "VQSR.mixedTest.input.vcf" +
                    " --output %s" +
                    " -tranchesFile " + getToolTestDataDir() + "VQSR.mixedTest.tranches" +
                    " -recalFile " + getToolTestDataDir() + "VQSR.mixedTest.recal.vcf",
                Arrays.asList(getToolTestDataDir() + "expected/applySNPAndIndelResult.vcf"));
        spec.executeTest("testApplyRecalibrationSnpAndIndelTogether", this);
    }

    @Test
    public void testApplyRecalibrationSnpAndIndelTogetherExcludeFiltered() throws Exception {
        ArgumentsBuilder args = new ArgumentsBuilder();
        File tempOut = createTempFile("testApplyRecalibrationSnpAndIndelTogetherExcludeFiltered", ".vcf");

        args.add("--variant");
        args.add(new File(getToolTestDataDir() + "VQSR.mixedTest.input.vcf"));
        args.add("-L");
        args.add("20:1000100-1000500");
        args.add("-mode");
        args.add("BOTH");
        args.add("--excludeFiltered");
        args.add("-ts_filter_level");
        args.add("90.0");
        args.add("-tranchesFile ");
        args.add(getToolTestDataDir() + "VQSR.mixedTest.tranches");
        args.add("-recalFile");
        args.add(getToolTestDataDir() + "VQSR.mixedTest.recal.vcf");
        args.addOutput(tempOut);

        runCommandLine(args);

        try (FeatureDataSource<VariantContext> featureSource = new FeatureDataSource<>(tempOut)) {
            for (VariantContext feature : featureSource) {
                // there should only be unfiltered records in the output VCF file
                Assert.assertTrue(feature.isNotFiltered());
            }
        }
    }

    @Test
    public void testApplyRecalibrationAlleleSpecificSNPmode() throws IOException {
        final String base =
                " -L 3:113005755-195507036" +
                " -mode SNP -AS" +
                " -ts_filter_level 99.7" +
                " --variant " + getToolTestDataDir() + "VQSR.AStest.input.vcf" +
                " --output %s" +
                " -tranchesFile " + getToolTestDataDir() + "VQSR.AStest.snps.tranches" +
                " -recalFile " + getToolTestDataDir() + "VQSR.AStest.snps.recal.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                base,
                Arrays.asList(getToolTestDataDir() + "expected/applySNPAlleleSpecificResult.vcf"));
        spec.executeTest("testApplyRecalibrationAlleleSpecificSNPmode", this);
    }

    @Test
    public void testApplyRecalibrationAlleleSpecificINDELmode() throws IOException {
        final String base =
                " -L 3:113005755-195507036" +
                " -mode INDEL -AS" +
                " -ts_filter_level 99.3" +
                " --variant " + getToolTestDataDir() + "VQSR.AStest.postSNPinput.vcf" +
                " --output %s" +
                " -tranchesFile " + getToolTestDataDir() + "VQSR.AStest.indels.tranches" +
                " -recalFile " + getToolTestDataDir() + "VQSR.AStest.indels.recal.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                base,
                Arrays.asList(getToolTestDataDir() + "expected/applyIndelAlleleSpecificResult.vcf"));
        spec.executeTest("testApplyRecalibrationAlleleSpecificINDELmode", this);
    }

}

