Create_data Proposal File Format
================================

The proposal file has an XML-like format that lays out the
relationship between a set of exposures.  The layout looks
like this:

.. literalinclude:: empty.prop

The file to be converted is put between the <base></base> tags,
and if a subarray is needed to be extracted from a full-frame exposure,
the name of the subarray can be put between the <subarray></subarray>
tags.  Finally, the type of exposure can be placed between the <exp_type>
</exp_type> tags.  The values of EXP_TYPE are:

+------------------+---------+----------------+---------+------------+
| MIRI             | NIRCAM  |  NIRSPEC       | NIRISS  |  FGS       |
+------------------+---------+----------------+---------+------------+
|MIR_IMAGE         |NRC_IMAGE|NRS_TASLIT      |NIS_IMAGE|FGS_IMAGE   |
+------------------+---------+----------------+---------+------------+
|MIR_TACQ          |NRC_TACQ |NRS_TACQ        |NIS_FOCUS|FGS_FOCUS   |
+------------------+---------+----------------+---------+------------+
|MIR_LYOT          |NRC_CORON|NRS_TACONFIRM   |NIS_DARK |FGS_SKYFLAT |
+------------------+---------+----------------+---------+------------+
|MIR_4QPM          |NRC_FOCUS|NRS_CONFIRM     |NIS_WFSS |FGS_INTFLAT |
+------------------+---------+----------------+---------+------------+
|MIR_LRS-FIXEDSLIT |NRC_DARK |NRS_FIXEDSLIT   |         |            |
+------------------+---------+----------------+---------+------------+
|MIR_LRS-SLITLESS  |NRC_FLAT |NRS_AUTOWAVECAL |         |            |  
+------------------+---------+----------------+---------+------------+
|MIR_MRS           |         |NRS_IFU         |         |            |
+------------------+---------+----------------+---------+------------+
|MIR_DARK          |         |NRS_MSA         |         |            |
+------------------+---------+----------------+---------+------------+
|MIR_FLAT          |         |NRS_AUTOFLAT    |         |            |
+------------------+---------+----------------+---------+------------+
|                  |         |NRS_DARK        |         |            |
+------------------+---------+----------------+---------+------------+
|                  |         |NRS_LAMP        |         |            |
+------------------+---------+----------------+---------+------------+

Sections of this file can be replicated to represent, for example,
all of the NIRCAM exposures from each of the 10 detectors at a single
pointing by just replicating the <detector></detector> blocks.
