/*
 * Copyright 2023 National Bank of Belgium
 * 
 * Licensed under the EUPL, Version 1.2 or – as soon they will be approved 
 * by the European Commission - subsequent versions of the EUPL (the "Licence");
 * You may not use this work except in compliance with the Licence.
 * You may obtain a copy of the Licence at:
 * 
 * https://joinup.ec.europa.eu/software/page/eupl
 * 
 * Unless required by applicable law or agreed to in writing, software 
 * distributed under the Licence is distributed on an "AS IS" basis,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the Licence for the specific language governing permissions and 
 * limitations under the Licence.
 */
package jdplus.dfm.base.core;

import internal.jdplus.dfm.base.core.CumulMeasurement;
import internal.jdplus.dfm.base.core.CumulatedVariationsMeasurement;
import internal.jdplus.dfm.base.core.LevelMeasurement;
import jdplus.dfm.base.api.MeasurementType;
import jdplus.toolkit.base.core.data.DataBlock;


/**
 * The IMeasurement interface represents the behaviour of a measurement equation
 * on a given factor (and its lags).
 *
 * @author Jean Palate
 */
public interface IDfmMeasurement {

    /**
     * The number of lags (nlags) implied by the measurement
     *
     * @return Lags in [t to t-nlags[ are used by the measurement in t
     */
    int getLength();

    /**
     * Fills the polynomial of the measurement (without the actual coefficient)
     * Typical values are: 1 [0 ... 0] 1 1 ... 1 1 2 3 2 1 [0 0... 0]
     *
     * @param z The buffer that will contain the polynomial. The length of the
     * buffer is defined by "getLength()".
     */
    void fill(DataBlock z);

    /**
     * Computes the product of the polynomial (z) with a block of data
     *
     * @param x The block of data. The length of the data is defined by
     * "getLength()".
     * @return returns z * x
     */
    double dot(DataBlock x);
    
       /**
     * Gets the default measurement for a given type.
     *
     * @param type The type of the measurement
     * @return For type C, returns (1...1) (length=12); for type CD, returns (1
     * 2 3 2 1); for type L, returns (1)
     */
    public static IDfmMeasurement measurement(final MeasurementType type) {
        return switch (type) {
            case YoY -> CumulMeasurement.MC12;
            case Q -> CumulatedVariationsMeasurement.MCD3;
            case M -> LevelMeasurement.ML;
            default -> null;
        };
    }

    public static MeasurementType getMeasurementType(final IDfmMeasurement m) {
        if (m instanceof CumulMeasurement) {
            return MeasurementType.YoY;
        } else if (m instanceof CumulatedVariationsMeasurement) {
            return MeasurementType.Q;
        } else if (m instanceof LevelMeasurement) {
            return MeasurementType.M;
        } else {
            return null;
        }
    }

}
