package org.broadinstitute.hellbender.tools.spark.sv;

import htsjdk.samtools.SAMUtils;
import htsjdk.samtools.fastq.FastqRecord;
import org.broadinstitute.hellbender.utils.read.GATKRead;
import org.testng.Assert;
import org.testng.annotations.Test;

/**
 * Created by shuang on 8/26/16.
 */
public class SVFastqUtilsUnitTest {

    @Test
    public void testFastqToGATKRead(){
        final FastqRecord customFastq = new FastqRecord("@HJYFJCCXX160204:4:1101:19004:15109/1 mapping=1:10275;151M",
                "CCCCCCAACCCCCAACCCTAACCCTAACCCCAACCCCAACCCCAACCCCAACCCTCACCCCAACCCCTAACCCCTAACCCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAACCCTAAC",
                "+",
                "7:99++/1?=**+/>@>:..0??9-.:?<8*./?>/*.0?=;*;;;:;*.8==;.)/=;2*.8==>;-.<>>>;.>0>><><80>>=>80>>><<0>>??=>>>?>====??>;>>>>=<==>>>>==>?>8>>??>>>>??>>>???>?=");

        final GATKRead read = SVFastqUtils.convertToRead(customFastq);

        Assert.assertEquals(read.getName(), customFastq.getReadHeader().split(" ")[0]);
        Assert.assertEquals(read.getBasesString(), customFastq.getReadString());
        Assert.assertEquals(SAMUtils.phredToFastq(read.getBaseQualities()), customFastq.getBaseQualityString());
        Assert.assertEquals(read.getAttributeAsString("RC"), "1");
    }
}
