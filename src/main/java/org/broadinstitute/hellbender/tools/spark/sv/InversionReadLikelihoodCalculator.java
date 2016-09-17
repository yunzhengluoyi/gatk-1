package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.variant.variantcontext.Allele;
import htsjdk.variant.variantcontext.VariantContext;
import org.apache.commons.lang3.StringUtils;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.genotyper.*;
import org.broadinstitute.hellbender.utils.pairhmm.PairHMM;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.Serializable;
import java.util.Arrays;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * RLC for inversion.
 */
final class InversionReadLikelihoodCalculator implements SVReadLikelihoodCalculator, Serializable{
    private static final long serialVersionUID = 1L;

    private ScoringScheme mScoring = new ScoringScheme();

    private PairHMM stolenHMM = PairHMM.Implementation.LOGLESS_CACHING.makeNewHMM();// FASTEST_AVAILABLE has dependency cannot be serialized....

    private int maxReadLength;
    private int maxAlleleLength;

    @Override
    public void close() {
        stolenHMM.close();
    }

    // TODO: right now doing nothing
    /**
     * {@inheritDoc}
     */
    @Override
    public List<GATKRead> preprocessReads(final List<GATKRead> reads){
        return SVReadLikelihoodCalculator.super.preprocessReads(reads);
    }

    /**
     * {@inheritDoc}
     */
    public ReadLikelihoods<SVDummyAllele> computeReadLikelihoods(final SampleList sampleList,
                                                                 final List<GATKRead> reads,
                                                                 final Map<String, List<GATKRead>> sample2Reads,
                                                                 final SVJunction junction){

        // get alleles, ref and alt, around the two breakpoints, from junction
        final List<SVDummyAllele> alleleList = junction.getAlleles();

        final ReadLikelihoods<SVDummyAllele> result = new ReadLikelihoods<>(sampleList, new IndexedAlleleList<>(alleleList), sample2Reads);
        if (reads.isEmpty()) {
            return result;
        }

        // (re)initialize the calculator only if necessary
        final int mrl = reads.stream().mapToInt(GATKRead::getLength).max().orElse(0);
        final int mal = junction.getAlleles().stream().mapToInt(Allele::length).max().orElse(0);
        if(!stolenHMM.isInitialized() || mrl > maxReadLength || mal > maxAlleleLength){
            maxReadLength = Math.max(mrl, maxReadLength);
            maxAlleleLength = Math.max(mal, maxAlleleLength);
            stolenHMM.initialize(maxReadLength, maxAlleleLength);
        }

        final String debugString = hackOfPairHMMComputeLog10Likelihoods(reads, result.sampleMatrix(0));
        if(GenotypingSingleSampleSVSpark.in_debug_state){
            final VariantContext vc = junction.getOriginalVC();
            junction.debugString = String.format("%s\n%s", junction.getOriginalVC().getID(), debugString);
        }

        return result;
    }

    // -----------------------------------------------------------------------------------------------
    // Shameless-copy from PairHMM class, for now
    // -----------------------------------------------------------------------------------------------

    private String hackOfPairHMMComputeLog10Likelihoods(final List<GATKRead> reads,
                                                        final LikelihoodMatrix<SVDummyAllele> logLikelihoods){
        // again, copied from HC
        final byte gcpDefault = (byte)10; // 10 is the default value for gcpHMM parameter in HC
        final Map<GATKRead, byte[]> gcp =  reads.stream().collect(Collectors.toMap(read -> read,read -> Utils.dupBytes(gcpDefault, read.getLength())));

        final List<SVDummyAllele> alleles = logLikelihoods.alleles();
        final int alleleCount = alleles.size();
        //debug array
        final double[] mLogLikelihoodArray = new double[reads.size() * alleleCount];

        int idx = 0;
        int readIndex = 0;
        for(final GATKRead read : reads){
            final byte[] readBases = read.getBases();
            final byte[] readQuals = read.getBaseQualities();
            final byte[] readInsQuals = ReadUtils.getBaseInsertionQualities(read);
            final byte[] readDelQuals = ReadUtils.getBaseDeletionQualities(read);
            final byte[] overallGCP = gcp.get(read);

            // peek at the next haplotype in the list (necessary to get nextHaplotypeBases, which is required for caching in the array implementation)
            final boolean isFirstHaplotype = true;
            for (int a = 0; a < alleleCount; a++) {
                final Allele allele = alleles.get(a);
                final byte[] alleleBases = allele.getBases();
                final byte[] nextAlleleBases = a == alleles.size() - 1 ? null : alleles.get(a + 1).getBases();
                final double lk = stolenHMM.computeReadLikelihoodGivenHaplotypeLog10(alleleBases, readBases, readQuals, readInsQuals, readDelQuals, overallGCP, isFirstHaplotype, nextAlleleBases);
                logLikelihoods.set(a, readIndex, lk);
                mLogLikelihoodArray[idx++] = lk;
            }
            readIndex++;
        }

        if(GenotypingSingleSampleSVSpark.in_debug_state){
            return formatDebugLikelihoodArray(mLogLikelihoodArray, alleles, reads.size());
        } else {
            return "";
        }
    }

    private static String formatDebugLikelihoodArray(final double[] readLL, final List<SVDummyAllele> alleles, final int rc){
        final StringBuilder builder = new StringBuilder();
        int i=0;
        for(final Allele a : alleles){
            builder.append(a.isReference() ? "*\t" : "ALT\t");
            final double[] l = Arrays.copyOfRange(readLL, i*rc, (i+1)*rc);
            for(int j=0; j<rc-1; ++j) builder.append(String.valueOf(l[j])).append(",");
            builder.append(String.valueOf(l[rc-1])).append("\n");
            ++i;
        }
        return builder.toString();
    }
}
