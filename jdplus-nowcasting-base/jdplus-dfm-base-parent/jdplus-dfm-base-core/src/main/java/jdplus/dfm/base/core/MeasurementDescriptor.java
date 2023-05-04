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

import static jdplus.dfm.base.core.IDfmMeasurement.getMeasurementType;
import jdplus.toolkit.base.api.data.DoubleSeq;



/**
 *
 * @author Jean Palate
 */
@lombok.Builder(builderClassName="Builder", toBuilder=true)
@lombok.Value
public class MeasurementDescriptor {

    public static final double C_DEF = .2;

    /**
     * Type of the measurement equation
     */
    private final IDfmMeasurement type;
    /**
     * Coefficients
     */
    private final DoubleSeq coefficient;
    /**
     * Variance of the measurement equation (>=0)
     */
    private double variance;

    public MeasurementDescriptor rescaleVariance(double c) {
        if (c == 1)
            return this;
        return new MeasurementDescriptor(type, coefficient, c*variance);
    }

    public MeasurementDescriptor withVariance(double var) {
        if (this.variance == var)
            return this;
        return new MeasurementDescriptor(type, coefficient, var);
    }

    public double getCoefficient(int pos) {
        return coefficient.get(pos);
    }

    @Override
    public String toString() {
        StringBuilder builder = new StringBuilder();
        for (int i = 0; i < coefficient.length(); ++i) {
            if (Double.isNaN(coefficient.get(i))) {
                builder.append('.');
            } else {
                builder.append(coefficient.get(i));
            }
            builder.append('\t');
        }
        builder.append(variance);
        return builder.toString();
    }

    public boolean isUsed(int fac) {
        return !Double.isNaN(coefficient.get(fac));
    }

    public boolean[] getUsedFactors() {
        boolean[] used = new boolean[coefficient.length()]; // DAVID: I THINK THIS LINE SHOULD BE COMMENTED (OTHERWISE IT CREATES A NEW USED)
        for (int i = 0; i < used.length; ++i) {
            used[i] = !Double.isNaN(coefficient.get(i));
        }
        return used;
    }

    public int getUsedFactorsCount() {
        
        return coefficient.count(x->!Double.isNaN(x));
    }

    MeasurementStructure getStructure() {
        return new MeasurementStructure(getMeasurementType(type), getUsedFactors());
    }
}
