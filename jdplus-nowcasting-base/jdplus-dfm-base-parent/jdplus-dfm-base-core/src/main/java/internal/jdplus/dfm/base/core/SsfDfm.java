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

import jdplus.dfm.base.core.MeasurementDescriptor;
import jdplus.dfm.base.core.TransitionDescriptor;
import jdplus.dfm.base.core.var.VarDescriptor;
import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.api.math.matrices.Matrix;
import jdplus.toolkit.base.core.ssf.multivariate.MultivariateSsf;
import org.checkerframework.checker.nullness.qual.NonNull;

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

    public MultivariateSsf unconditionalSsf(TransitionDescriptor vdesc, MeasurementDescriptor[] mdesc, int extendedLags) {
        int nlx = extendedLags(Math.max(extendedLags, vdesc.getNlags()), mdesc);
        int nf = vdesc.getNfactors();
        Dynamics dyn = Dynamics.of(vdesc.getCoefficients(), vdesc.getInnovationsVariance(), nlx);
        Initialization initialization = Initialization.unconditional(dyn);
        Measurements m = Measurements.of(nf, nlx, mdesc);
        return new MultivariateSsf(initialization, dyn, m);
    }

    public MultivariateSsf userSsf(TransitionDescriptor vdesc, MeasurementDescriptor[] mdesc, DoubleSeq a0, @NonNull Matrix V0) {
        int dim = V0.getRowsCount();
        Dynamics dyn = Dynamics.of(vdesc.getCoefficients(), vdesc.getInnovationsVariance(), dim);
        Initialization initialization = Initialization.user(a0, V0);
        Measurements m = Measurements.of(vdesc.getNfactors(), dim, mdesc);
        return new MultivariateSsf(initialization, dyn, m);
    }

    public MultivariateSsf zeroSsf(TransitionDescriptor vdesc, MeasurementDescriptor[] mdesc, int extendedLags) {
        int nlx = extendedLags(Math.max(extendedLags, vdesc.getNlags()), mdesc);
        int nf = vdesc.getNfactors();
        Dynamics dyn = Dynamics.of(vdesc.getCoefficients(), vdesc.getInnovationsVariance(), nlx);
        Initialization initialization = Initialization.zero(nlx);
        Measurements m = Measurements.of(nf, nlx, mdesc);
        return new MultivariateSsf(initialization, dyn, m);
    }

}
