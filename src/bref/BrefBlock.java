/*
 * Copyright 2021-2023 Brian L. Browning
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
package bref;

import beagleutil.ChromIds;
import blbutil.Const;

/**
 * <p>Class {@code BrefBlock} represents starting chromosome coordinates and
 * file offset for the start of a binary reference format (bref) data block.
 * </p>
 *
 * <p>Instances of class {@code BrefBlock} are immutable.
 * </p>
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public class BrefBlock {

    private final int chromIndex;
    private final int pos;
    private final long offset;

    /**
     * Constructs a {@code BrefBlock} for the specified data.  It is the
     * caller's responsibility to ensure the consistency of the
     * constructor parameters.
     *
     * @param chromIndex the chromosome index
     * @param pos the starting chromosome position
     * @param offset the file offset in bytes for the bref data block
     */
    public BrefBlock(int chromIndex, int pos, long offset) {
        this.chromIndex = chromIndex;
        this.pos = pos;
        this.offset = offset;
    }

    /**
     * Returns the chromosome index of the first marker in this bref block.
     * @return the chromosome index of the first marker in this bref block
     */
    public int chromIndex() {
        return chromIndex;
    }

    /**
     * Returns the chromosome position of the first marker in this bref block.
     * @return the chromosome position of the first marker in this bref block
     */
    public int pos() {
        return pos;
    }

    /**
     * Returns the file offset of the first marker in this bref block.
     * @return the file offset of the first marker in this bref block
     */
    public long offset() {
        return offset;
    }

    /**
     * Returns a string description of {@code this}.  The exact details
     * of the representation are unspecified and subject to change.
     * @return  a string description of {@code this}
     */
    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder(30);
        sb.append('[');
        sb.append(ChromIds.instance().id(chromIndex));
        sb.append(Const.tab);
        sb.append(pos);
        sb.append(Const.tab);
        sb.append(offset);
        sb.append(']');
        return sb.toString();
    }
}
