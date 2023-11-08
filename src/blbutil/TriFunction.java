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
package blbutil;

/**
 * <p>Class {@code TriFunction} represents a function that accepts three
 * arguments and produces a result.  This functional interface that
 * is analogous to the {@code java.util.function.BiFunction} interface.
 *
 * @param <T> the type of the first argument to the function
 * @param <U> the type of the second argument to the function
 * @param <V> the type of the third argument to the function
 * @param <R> the type of the result of the function
 *
 * @author Brian L. Browning {@code <browning@uw.edu>}
 */
public interface TriFunction<T, U, V, R> {

    /**
     * Applies this function to the specified arguments.
     *
     * @param t the first function argument
     * @param u the second function argument
     * @param v the third function argument
     * @return the function result
     * @throws NullPointerException if any function argument is {@code null}
     */

    R apply(T t, U u, V v);
}
