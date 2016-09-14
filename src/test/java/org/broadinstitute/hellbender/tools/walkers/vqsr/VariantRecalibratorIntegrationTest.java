package org.broadinstitute.hellbender.tools.walkers.vqsr;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.testng.Assert;
import org.testng.annotations.BeforeMethod;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.stream.Stream;

public class VariantRecalibratorIntegrationTest extends CommandLineProgramTest {

    @Override
    public String getTestedClassName() {
        return VariantRecalibrator.class.getSimpleName();
    }

    public String getToolTestDataDir(){
        return publicTestDir + "org/broadinstitute/hellbender/tools/walkers/VQSR/";
    }

    public String getLargeVQSRTestDataDir(){
        return largeFileTestDir + "VQSR/";
    }

    @BeforeMethod
    public void initializeVariantRecalTests() {
        //Reset the RNG in order to align it with the initial state of the GATK3 RNG at the time
        //the tests start running, which we want to do in order to get the same results produced by
        //GATK3. Note that this means the results of running these tests will be different if they
        //are run manually outside of the test framework.
        logger.info("Initializing VQSR tests/resetting random number generator");
        Utils.resetRandomGenerator();
    }

    @DataProvider(name="VarRecalSNP")
    public Object[][] getVarRecalSNPData() {
        return new Object[][]{
            {
                new String[] {
                    "--variant",
                    getLargeVQSRTestDataDir() + "phase1.projectConsensus.chr20.1M-10M.raw.snps.vcf",
                    "-L","20:1,000,000-10,000,000",
                    "--resource",
                    "known,known=true,prior=10.0:" + getLargeVQSRTestDataDir() + "dbsnp_132_b37.leftAligned.20.1M-10M.vcf",
                    "--resource",
                    "truth_training1,truth=true,training=true,prior=15.0:" + getLargeVQSRTestDataDir() + "sites_r27_nr.b37_fwd.20.1M-10M.vcf",
                    "--resource",
                    "truth_training2,training=true,truth=true,prior=12.0:" + getLargeVQSRTestDataDir() + "Omni25_sites_1525_samples.b37.20.1M-10M.vcf",
                    "-an", "QD", "-an", "HaplotypeScore", "-an", "HRun",
                    "--trustAllPolymorphic", // for speed
                    "-mode", "SNP"
                }
            },
        };
    }

    @Test(dataProvider = "VarRecalSNP")
    public void testVariantRecalibratorSNP(final String[] params) throws IOException {
        //NOTE: The number of iterations required to ensure we have enough negative training data to proceed,
        //as well as the test results themselves, are both very sensitive to the state of the random number
        //generator at the time the tool starts to execute. Sampling a single integer from the RNG at the
        //start aligns the initial state of the random number generator with the initial state of the GATK3
        //random number generator at the time the tests start executing (both RNGs use the same seed, but
        //the pathway to the test in GATK3 results in an additional integer being sampled), thus allowing the
        //results to be easily compared against those produced by GATK3. This also happens to allow this
        //test to succeed on the first iteration (max_attempts=1) which we want for test performance reasons.
        //Failing to call nextInt here would result in the model failing on the first 3 attempts (it would
        //succeed on the 4th), and the results would not match those produced by GATK3 on these same inputs.
        //
        //Also note that due to this RNG conditioning, this test will produce different results when manually
        //outside of the test framework.
        @SuppressWarnings("unused")
        final int hack = Utils.getRandomGenerator().nextInt();

        // use an ArrayList - ArgumentBuilder tokenizes using the "=" in the resource args
        List<String> args = new ArrayList<>(params.length);
        Stream.of(params).forEach(arg -> args.add(arg));

        File recalOut = createTempFile("testVarRecalSnp", ".vcf");
        File tranchesOut = createTempFile("testVarRecalSnp", ".txt");
        args.addAll(addTempFileArgs(recalOut, tranchesOut));

        final VariantRecalibrator varRecalTool = new VariantRecalibrator();
        Assert.assertEquals(varRecalTool.instanceMain(args.toArray(new String[args.size()])), true);

        // the expected vcf is not in the expected dir because its used
        // as input for the ApplyVQSR test
        IntegrationTestSpec.assertEqualTextFiles(recalOut, new File(getLargeVQSRTestDataDir() + "snpRecal.vcf"));
        IntegrationTestSpec.assertEqualTextFiles(tranchesOut, new File(getLargeVQSRTestDataDir() + "expected/SNPTranches.txt"));
    }

    @Test(dataProvider = "VarRecalSNP")
    public void testVariantRecalibratorSNPMaxAttempts(final String[] params) throws IOException {
        // For this test, we deliberately *DON'T* sample a single random int as above; this causes
        // the tool to require 4 attempts to acquire enough negative training data to succeed

        // use an ArrayList - ArgumentBuilder tokenizes using the "=" in the resource args
        List<String> args = new ArrayList<>(params.length);
        Stream.of(params).forEach(arg -> args.add(arg));
        File recalOut = createTempFile("testVarRecalMaxAttempts", ".vcf");
        File tranchesOut = createTempFile("testVarRecalMaxAttempts", ".txt");
        args.addAll(addTempFileArgs(recalOut, tranchesOut));

        args.add("--max_attempts");
        args.add("4"); // it takes for for this test to wind up with enough training data

        final VariantRecalibrator varRecalTool = new VariantRecalibrator();
        Assert.assertEquals(varRecalTool.instanceMain(args.toArray(new String[args.size()])), true);
        Assert.assertEquals(varRecalTool.max_attempts, 4);
    }

    private List<String> addTempFileArgs(final File recalOutFile, final File tranchesOutFile) {
        List<java.lang.String> args = new ArrayList<>(2);
        args.add("--output");
        args.add(recalOutFile.getAbsolutePath());
        args.add("--tranches_file");
        args.add(tranchesOutFile.getAbsolutePath());
        return args;
    }

    @Test
    public void testVariantRecalibratorIndel() throws IOException {
        @SuppressWarnings("unused")
        final int hack = Utils.getRandomGenerator().nextInt();
        final String inputFile = getLargeVQSRTestDataDir() + "combined.phase1.chr20.raw.indels.filtered.sites.1M-10M.vcf";

        final IntegrationTestSpec spec = new IntegrationTestSpec(
                " --resource known,known=true,prior=10.0:" + getLargeVQSRTestDataDir() + "dbsnp_132_b37.leftAligned.20.1M-10M.vcf" +
                " --resource truth_training,training=true,truth=true,prior=15.0:" + getLargeVQSRTestDataDir() + "ALL.wgs.indels_mills_devine_hg19_leftAligned_collapsed_double_hit.sites.20.1M-10M.vcf" +
                " --variant " + inputFile +
                " -L 20:1,000,000-10,000,000" +
                " -an QD -an ReadPosRankSum -an HaplotypeScore" +
                " -mode INDEL -mG 3" +
                " --trustAllPolymorphic" + // for speed
                " --output %s" +
                " -tranchesFile %s",
                Arrays.asList(
                        // the "expected" vcf is not in the expected dir because its used
                        // as input for the ApplyVQSR test
                        getLargeVQSRTestDataDir() + "indelRecal.vcf",
                        getLargeVQSRTestDataDir() + "expected/indelTranches.txt"));
        spec.executeTest("testVariantRecalibratorIndel"+  inputFile, this);
    }

}

