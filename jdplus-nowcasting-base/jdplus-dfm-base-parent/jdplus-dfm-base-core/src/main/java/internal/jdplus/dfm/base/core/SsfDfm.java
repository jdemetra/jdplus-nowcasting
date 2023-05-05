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
package internal.jdplus.dfm.base.core;

import jdplus.dfm.base.core.DynamicFactorModel;
import jdplus.dfm.base.core.MeasurementDescriptor;
import jdplus.dfm.base.core.var.VarDescriptor;
import jdplus.toolkit.base.api.math.matrices.Matrix;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;
import jdplus.toolkit.base.core.ssf.multivariate.MultivariateSsf;

/**
 *
 * @author Jean Palate
 */
@lombok.experimental.UtilityClass
public class SsfDfm {

    private int extendedLags(int nlags, MeasurementDescriptor[] mdesc) {
        int nlx = nlags;
        for (int i = 0; i < mdesc.length; ++i) {
            int n = mdesc[i].getType().getLength();
            if (nlx < n) {
                nlx = n;
            }
        }
        return nlx;
    }

    private static FastMatrix convert(Matrix M) {
        int nf = M.getRowsCount(), nc = M.getColumnsCount();
        int nl = nc / nf;
        FastMatrix C = FastMatrix.make(nf, nf * nl);
        // M columns arrange by lags, C columns arranged by factors
        for (int l = 0, k = 0; l < nl; ++l) {
            for (int f = 0; f < nf; ++f, ++k) {
                C.column(l + f * nl).copy(M.column(k));
            }
        }
        return C;
    }

    public MultivariateSsf of(DynamicFactorModel dfm, int extendedLags) {
        VarDescriptor var = dfm.getVar();
        MeasurementDescriptor[] mdesc = dfm.getMeasurements().toArray(MeasurementDescriptor[]::new);
        int nxlags = extendedLags(Math.max(extendedLags, var.getNlags()), mdesc);
        int nf = var.getNfactors();
        Dynamics dyn = Dynamics.of(convert(var.getCoefficients()), FastMatrix.of(var.getInnovationsVariance()), nxlags);
        Initialization initialization = switch (var.getInitialization()) {
            case Zero ->
                Initialization.zero(var.getInnovationsVariance(), nxlags);
            case Unconditional ->
                Initialization.unconditional(dyn);
            default ->
                throw new IllegalArgumentException();
        };
        Measurements m = Measurements.of(nf, nxlags, mdesc);
        return new MultivariateSsf(initialization, dyn, m);
    }

}
