package org.broadinstitute.hellbender.tools.spark.sv;

import com.google.common.annotations.VisibleForTesting;
import htsjdk.variant.variantcontext.Genotype;
import htsjdk.variant.variantcontext.GenotypeBuilder;
import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.variantcontext.VariantContextBuilder;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.genotyper.AlleleList;
import org.broadinstitute.hellbender.utils.variant.GATKVariantContextUtils;
import scala.Tuple3;

import java.io.IOException;
import java.io.Serializable;
import java.util.*;
import java.util.stream.Collectors;

/**
 * An {@code SVJunction} represents a novel adjacency (or in the case of insertion, disjoint)
 * between TWO (or in the case of insertion, ONE) genomic locations on the reference created by a genomic structural event.
 */
abstract class SVJunction implements Serializable{
    private static final long serialVersionUID = 1L;

    // -----------------------------------------------------------------------------------------------
    // Common base fields for all types of SV's. All should be set by parsing caller output VCF.
    // -----------------------------------------------------------------------------------------------

    /**
     * The variant context that was constructed by caller.
     * Will not be modified but may be used for constructing new VC's based on it.
     */
    protected final VariantContext vc;

    protected final GATKSVVCFHeaderLines.SVTYPES type;

    /**
     * 5'- and 3'- breakpoints locations. Both will be a single base location.
     * For insertion these two would be the same.
     */
    protected final SimpleInterval fiveEndBPLoc;
    protected final SimpleInterval threeEndBPLoc;

    protected byte whichEnd; // only 4 possible states 0: neither (why called?); 1: 5'-end; 2: 3'-end; 3: both (error?)

    /**
     * A list of assembly contigs that caller thinks support this junction.
     */
    protected final List<LocalAssemblyContig> contigsInformation;

    /**
     * Number of high quality mappings (criteria decided by caller).
     * Less than or equal to the length of {@link #contigsInformation}.
     */
    protected final int numHQMappings;

    /**
     * Maximal alignment length out of the contigs that support this breakpoint
     */
    protected final int maxAlignLength;

    /**
     * Length of the SV event w.r.t the reference as determined by caller.
     * Using type int because an sv length greater than int.max (2 billion) should have been filtered out by caller.
     */
    protected final int svLength;

    /**
     * (Micro-)Homology around supposedly two breakpoints of the same SV event.
     * null when the VC doesn't contain this attribute.
     */
    protected final byte[] maybeNullHomology;

    /**
     * Inserted sequence around an structural variation breakpoint.
     * null when the VC doesn't contain this attribute.
     */
    protected final byte[] maybeNullInsertedSeq;


    String debugString; // a place to put debugging information

    // -----------------------------------------------------------------------------------------------
    // Fields in this block must be computed, not parsed from input VC record
    // -----------------------------------------------------------------------------------------------

    // 5'- and 3'- end reference intervals that the reference alleles--one for each breakpoint--reside in
    protected Tuple3<String, String, int[]> referenceWindows = new Tuple3<>("", "", new int[] {-1, -1, -1, -1});

    /**
     * List of alleles for this junction.
     * Actual instances must return an even number for {@link AlleleList#numberOfAlleles()},
     * because of two breakpoints (5'- and 3'- end), except for insertion where the two breakpoints coincide.
     */
    protected List<SVDummyAllele> alleleList = new ArrayList<>();


    private int[] plVec = new int[] {0, 0, 0};

    // -----------------------------------------------------------------------------------------------
    // Common interfaces (mostly trivia)
    // -----------------------------------------------------------------------------------------------

    final VariantContext getOriginalVC(){ return vc;}

    /**
     * Get ids of putative breakpoints identified by {@link FindBreakpointEvidenceSpark}, that the caller believes supports this SV junction.
     */
    final List<Long> getAssemblyIDs(){
        return contigsInformation.stream().map(info -> info.assemblyID).collect(Collectors.toList());
    }

    final List<String> getContigIDs(){
        return contigsInformation.stream().map(info -> info.contigID).collect(Collectors.toList());
    }

    /**
     * Get the positions on reference this SV junction brought together (or set apart if it is an insertion).
     */
    final List<SimpleInterval> getLeftAlignedBreakpointLocations(){
        return Arrays.asList(fiveEndBPLoc, threeEndBPLoc);
    }

    /**
     * @return alleles around breakpoints of this junction. The reference alleles are guaranteed to come first.
     */
    final List<SVDummyAllele> getAlleles(){
        return alleleList;
    }

    final Tuple3<String, String, int[]> getReferenceWindows(){
        return referenceWindows;
    }

    /**
     * Since GL is most likely deprecated, no GL accessor is provided, only PL.
     */
    final int[] getPL(){
        return Arrays.copyOfRange(plVec, 0, plVec.length);
    }

    final SVJunction setPL(final int[] pls) {
        plVec = Arrays.copyOfRange(pls, 0, pls.length);
        return this;
    }

    /**
     * Add new attributes to the VC that was used for constructing this junction record
     * to the VC stored in the junction that just went under genotyping and return it.
     *
     * Note that FT will not be handled here.
     * TODO: works only for single sample
     */
    final VariantContext getGenotypedVC(){

        final Genotype tobeDetermined = new GenotypeBuilder(GenotypingSingleSampleSVSpark.testSampleName, GATKVariantContextUtils.noCallAlleles(2))
                .PL(getPL())
                .noDP()
                .noAD()
                .phased(false)
                .make();

        return new VariantContextBuilder(this.vc).genotypes(GenotypingSingleSampleSVSpark.inferGenotypeFromPL(tobeDetermined, vc.getAlleles())).make();
    }

    // -----------------------------------------------------------------------------------------------
    // Abstracts
    // -----------------------------------------------------------------------------------------------

    // TODO: figure out why using the list contigsInformation would lead to two same junctions being evaluated to be
    //       different, one plain one constructed by Spark
    @Override
    public boolean equals(final Object o) {
        if (this == o) return true;
        if (o == null || getClass() != o.getClass()) return false;

        final SVJunction junction = (SVJunction) o;

        if (whichEnd != junction.whichEnd) return false;
        if (svLength != junction.svLength) return false;
        if (type != junction.type) return false;
        if (!fiveEndBPLoc.equals(junction.fiveEndBPLoc)) return false;
        if (!threeEndBPLoc.equals(junction.threeEndBPLoc)) return false;
        if (!vc.getAttributeAsString(GATKSVVCFHeaderLines.ASSEMBLY_IDS, " ").equals(junction.vc.getAttributeAsString(GATKSVVCFHeaderLines.ASSEMBLY_IDS, " "))) return false;
        return vc.getAttributeAsString(GATKSVVCFHeaderLines.CONTIG_IDS, " ").equals(junction.vc.getAttributeAsString(GATKSVVCFHeaderLines.CONTIG_IDS, " "));
    }

    @Override
    public int hashCode() {
        int result = type.ordinal();
        result = 31 * result + fiveEndBPLoc.hashCode();
        result = 31 * result + threeEndBPLoc.hashCode();
        result = 31 * result + (int) whichEnd;
        result = 31 * result + svLength;
        result = 31 * result + vc.getAttributeAsString(GATKSVVCFHeaderLines.ASSEMBLY_IDS, " ").hashCode();
        result = 31 * result + vc.getAttributeAsString(GATKSVVCFHeaderLines.CONTIG_IDS, " ").hashCode();
        return result;
    }

    /**
     * Given a {@link VariantContext} (most likely from a parse of VCF spewed by the caller, i.e. DISCOVERY mode),
     * the ctor parses the VC and initializes (or input) appropriate fields.
     * @param vc
     */
    protected SVJunction(final VariantContext vc,
                         final Map<Long, List<LocalAssemblyContig>> assemblyID2ContigsAlignRecords){

        // simply copy from caller output VCF
        this.vc = vc;
        // TODO: is this really necessary?
        type = Enum.valueOf(GATKSVVCFHeaderLines.SVTYPES.class, (String)vc.getAttribute(GATKSVVCFHeaderLines.SVTYPE));

        fiveEndBPLoc = new SimpleInterval( vc.getContig() , vc.getStart() , vc.getStart() );
        threeEndBPLoc = new SimpleInterval( vc.getContig() , vc.getEnd() , vc.getEnd() );

        contigsInformation = new ArrayList<>(400); // blunt guess
        fillContigsInformation(contigsInformation, assemblyID2ContigsAlignRecords);

        numHQMappings = vc.getAttributeAsInt(GATKSVVCFHeaderLines.HQ_MAPPINGS, -1);
        maxAlignLength = vc.getAttributeAsInt(GATKSVVCFHeaderLines.MAX_ALIGN_LENGTH, -1);
        svLength = vc.getAttributeAsInt(GATKSVVCFHeaderLines.SVLEN, -1);
        if(numHQMappings==-1 || maxAlignLength==-1 || svLength==-1){
            throw new IllegalArgumentException("Seemingly broken VCF with any of "  +
                    GATKSVVCFHeaderLines.HQ_MAPPINGS + GATKSVVCFHeaderLines.MAX_ALIGN_LENGTH + GATKSVVCFHeaderLines.SVLEN
                    + " missing at site: " + vc.getContig() + ":" + vc.getStart());
        }

        final Object maybeNullHomologyObj = vc.getAttribute(GATKSVVCFHeaderLines.HOMOLOGY);
        maybeNullHomology = maybeNullHomologyObj ==null? null : ((String) maybeNullHomologyObj).getBytes();
        final Object maybeNullInsertedSeqObj = vc.getAttribute(GATKSVVCFHeaderLines.INSERTED_SEQUENCE);
        maybeNullInsertedSeq = maybeNullInsertedSeqObj ==null? null : ((String) maybeNullInsertedSeqObj).getBytes();
        // lazy work is over, no more copying, compute
    }

    protected abstract void setWhichEnd();

    /**
     * Construct windows around the two breakpoints identified.
     * @return the returned struct has string representations of the two reference contig names for the two breakpoints,
     *         first 5'-End then 3'-End, the int array is a length-4 array indicating the starting and ending points of the
     *         two windows, all inclusive.
     */
    protected abstract Tuple3<String, String, int[]> constructReferenceWindows();

    protected void constructAlleles(final ReferenceMultiSource reference){
        referenceWindows = constructReferenceWindows();

        final List<SVDummyAllele> alleles = constructReferenceAlleles(reference);
        alleles.addAll(constructAlternateAlleles(reference));
        alleleList = alleles;
    }

    /**
     * Construct reference alleles around the two breakpoints of the identified sv junction.
     *
     * By default, it simply extract bases in the windows marked by {@link #constructReferenceWindows()}.
     *
     * @return a list of reference alleles, with 5'-end being the first and 3'-end being the second. In case of inversion, only one will be .
     */
    @VisibleForTesting
    protected List<SVDummyAllele> constructReferenceAlleles(final ReferenceMultiSource reference){
        try{
            final SimpleInterval fiveEndWindow = new SimpleInterval(referenceWindows._1(), referenceWindows._3()[0], referenceWindows._3()[1]);
            final SimpleInterval threeEndWindow = new SimpleInterval(referenceWindows._2(), referenceWindows._3()[2], referenceWindows._3()[3]);

            return new ArrayList<>(Arrays.asList(new SVDummyAllele(reference.getReferenceBases(null, fiveEndWindow).getBases(), true),
                                 new SVDummyAllele(reference.getReferenceBases(null, threeEndWindow).getBases(), true)));

        } catch (final IOException ioex){
            throw new GATKException("Cannot resolve reference for constructing ref allele for junction.");
        }
    }

    protected abstract List<SVDummyAllele> constructAlternateAlleles(final ReferenceMultiSource reference);

    // -----------------------------------------------------------------------------------------------
    // Utilities
    // -----------------------------------------------------------------------------------------------

    /**
     * Extract from the input map {@code assemblyID2ContigsAlignRecords} informative contigs for this particular junction,
     * and use the extracted contigs to fill {@code tobeFilled}.
     */
    private void fillContigsInformation(final List<LocalAssemblyContig> tobeFilled,
                                        final Map<Long, List<LocalAssemblyContig>> assemblyID2ContigsAlignRecords){

        final Set<Long> asmids = Arrays.stream( vc.getAttributeAsString(GATKSVVCFHeaderLines.ASSEMBLY_IDS, "").replace("[", "").replace("]", "").split(", ") ).map(Long::valueOf).collect(Collectors.toSet());
        final Set<String> cids = Arrays.stream( vc.getAttributeAsString(GATKSVVCFHeaderLines.CONTIG_IDS, "").replace("[", "").replace("]", "").split(", ") ).map(s -> s.replace(">", "")).collect(Collectors.toSet());
        for(final Long asmID : asmids){
            tobeFilled.addAll(assemblyID2ContigsAlignRecords.get(asmID).stream().filter( contig -> cids.contains(contig.contigID)).collect(Collectors.toList()));
        }
    }
}
