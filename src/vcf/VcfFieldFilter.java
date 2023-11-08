/*
 * Copyright 2021 Brian L. Browning
 *
 * This file is part of the flare program.
 *
 * Licensed under the Apache License, Version 2.0 (the License);
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an AS IS BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
package vcf;

import blbutil.Const;
import blbutil.StringUtil;

/**
 * <p>Class {@code VcfFieldFilter} is a filter for a VCF record's ID, QUAL, 
 * FILTER, and INFO subfields.</p>
 *
 * <p>Instances of class {@code VcfFieldFilter} are immutable.</p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public final class VcfFieldFilter {

    private final boolean storeId;
    private final boolean storeQual;
    private final boolean storeFilter;
    private final boolean storeInfo;

    /**
     * Constructs a new {@code VcfFieldFilter} instance from the specified
     * data.  The INFO/END subfield always passes this filter.  Thus if
     * {@code (storeInfo == false)} and a VCF record INFO/END field is present,
     * the INFO/END field will pass this filter.
     * @param storeId {@code true} if a non-missing VCF ID field will be stored
     * @param storeQual {@code true} if non-missing VCF QUAL field will be stored
     * @param storeFilter {@code true} if a non-missing VCF FILTER field will be stored
     * @param storeInfo {@code true} if a non-missing VCF INFO field will be stored
     */
    public VcfFieldFilter(boolean storeId, boolean storeQual,
            boolean storeFilter, boolean storeInfo) {
        this.storeId = storeId;
        this.storeQual = storeQual;
        this.storeFilter = storeFilter;
        this.storeInfo = storeInfo;
    }


    /**
     * Appends the VCF ID subfields that pass this filter to the
     * specified {@code StringBuilder}.  Returns {@code true} if one or
     * more ID subfields are appended to the specified
     * {@code StringBuilder}, and returns {@code false} otherwise.
     * @param vcfRec a VCF record
     * @param start the index (inclusive) of the first character in the VCF
     * ID field
     * @param end the index (exclusive) of the last character in the VCF ID
     * field
     * @param sb the {@code StringBuilder} to which a non-missing ID field
     * may be appended
     * @return {@code true} if one or more ID subfields were appended
     * to the specified {@code StringBuilder}
     *
     * @throws IllegalArgumentException if an ID subfield passes this filter,
     * and {@code (start < 0)}, or {@code (start > end)},
     * or {@code (end > vcfRec.length())}
     * @throws NullPointerException if {@code (vcfRec == null) || (sb == null)}
     */
    public boolean storeId(String vcfRec, int start, int end, StringBuilder sb) {
        if (sb == null) {
            throw new IllegalArgumentException(StringBuilder.class.toString());
        }
        if (storeId) {
            if ((end-start)==1 && vcfRec.charAt(start)==Const.MISSING_DATA_CHAR) {
                return false;
            }
            else {
                sb.append(vcfRec, start, end);
                return true;
            }
        }
        else {
            return false;
        }
    }

    /**
     * Appends the VCF QUAL field to the specified {@code StringBuilder} if
     * the QUAL field passes this filter.
     * Returns {@code true} if the QUAL field is appended to the
     * specified {@code StringBuilder}, and returns {@code false} otherwise.
     * @param vcfRec a VCF record
     * @param start the index (inclusive) of the first character in the VCF QUAL
     * field
     * @param end the index (exclusive) of the last character in the VCF QUAL
     * field
     * @param prependTab {@code true} if an appended field should be
     * preceded by a tab ('\t')
     * @param sb the {@code StringBuilder} to which the QUAL field may be
     * appended
     * @return {@code true} if the QUAL field was appended to the specified
     * {@code StringBuilder}
     *
     * @throws IllegalArgumentException if the QUAL field passes this filter
     * and {@code (start < 0)}, or {@code (start > end)},
     * or {@code (end > vcfRec.length())}
     * @throws NullPointerException if {@code (vcfRec == null) || (sb == null)}
     */
    public boolean storeQual(String vcfRec, int start, int end,
            boolean prependTab, StringBuilder sb) {
        if (sb == null) {
            throw new IllegalArgumentException(StringBuilder.class.toString());
        }
        if (storeQual) {
            if ((end-start)==1 && vcfRec.charAt(start)==Const.MISSING_DATA_CHAR) {
                return false;
            }
            else {
                if (prependTab) {
                    sb.append(Const.tab);
                }
                sb.append(vcfRec, start, end);
                return true;
            }
        }
        else {
            return false;
        }
    }

    /**
     * Appends the VCF FILTER subfields that pass this filter to the
     * specified {@code StringBuilder}.  Returns {@code true} if one or
     * more FILTER subfields are appended to the specified
     * {@code StringBuilder}, and returns {@code false} otherwise.
     * @param vcfRec a VCF record
     * @param start the index (inclusive) of the first character in the VCF
     * FILTER field
     * @param end the index (exclusive) of the last character in the VCF FILTER
     * field
     * @param prependTab {@code true} if appended subfields should be
     * preceded by a tab ('\t')
     * @param sb the {@code StringBuilder} to which FILTER subfields may be
     * appended
     * @return {@code true} if one or more FILTER subfields were appended
     * to the specified {@code StringBuilder}
     *
     * @throws IllegalArgumentException if a FILTER subfield passes this filter,
     * and {@code (start < 0)}, or {@code (start > end)},
     * or {@code (end > vcfRec.length())}
     * @throws NullPointerException if {@code (vcfRec == null) || (sb == null)}
     */
    public boolean storeFilter(String vcfRec, int start, int end,
            boolean prependTab, StringBuilder sb) {
        if (sb == null) {
            throw new IllegalArgumentException(StringBuilder.class.toString());
        }
        if (storeFilter) {
            if ((end-start)==1 && vcfRec.charAt(start)==Const.MISSING_DATA_CHAR) {
                return false;
            }
            else {
                if (prependTab) {
                    sb.append(Const.tab);
                }
                sb.append(vcfRec, start, end);
                return true;
            }
        }
        else {
            return false;
        }
    }

    /**
     * Appends the VCF INFO subfields that pass this filter to the
     * specified {@code StringBuilder}.Returns {@code true} if one or
     * more INFO subfields are appended to the specified
     * {@code StringBuilder}, and returns {@code false} otherwise.
     * @param vcfRec a VCF record
     * @param start the index (inclusive) of the first character in the VCF
     * INFO field
     * @param end the index (exclusive) of the last character in the VCF INFO
     * field
     * @param prependTab {@code true} if appended subfields should be
     * preceded by a tab ('\t')
     * @param sb the {@code StringBuilder} to which INFO subfields may be
     * appended
     * @return {@code true} if one or more INFO subfields were appended
     * to the specified {@code StringBuilder}
     *
     * @throws IllegalArgumentException if an INFO subfield passes this filter,
     * and {@code (start < 0)}, or {@code (start > end)},
     * or {@code (end > vcfRec.length())}
     * @throws NullPointerException if {@code (vcfRec == null) || (sb == null)}
     */
    public boolean storeInfo(String vcfRec, int start, int end,
            boolean prependTab, StringBuilder sb) {
        if (sb == null) {
            throw new IllegalArgumentException(StringBuilder.class.toString());
        }
        if (storeInfo) {
            if ((end-start)==1 && vcfRec.charAt(start)==Const.MISSING_DATA_CHAR) {
                return false;
            }
            else {
                if (prependTab) {
                    sb.append(Const.tab);
                }
                sb.append(vcfRec, start, end);
                return true;
            }
        }
        else {
            String endSubfield = endSubfield(vcfRec, start, end);
            if (endSubfield == null) {
                return false;
            }
            else {
                if (prependTab) {
                    sb.append(Const.tab);   // tab before INFO field
                }
                sb.append(endSubfield);
                return true;
            }
        }
    }

    /*
     * Returns the first INFO/END subfield, or {@code null} if there
     * is no INFO/END subfield.
     */
    private static String endSubfield(String vcfRec, int infoStart, int infoEnd) {
        String infoField = vcfRec.substring(infoStart, infoEnd);
        String[] fields = StringUtil.getFields(infoField, Const.semicolon);
        for (String field : fields) {
            if (field.startsWith("END=")) {
                return field;
            }
        }
        return null;
    }

    /**
     * Returns a string description of {@code this}.  The exact details of
     * the description of are unspecified and subject to change.
     *
     * @return a string description of {@code this}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(50);
        sb.append("ID=");
        sb.append(storeId);
        sb.append(" QUAL=");
        sb.append(storeQual);
        sb.append(" FILTER=");
        sb.append(storeFilter);
        sb.append(" INFO=");
        sb.append(storeInfo);
        return sb.toString();
    }
}