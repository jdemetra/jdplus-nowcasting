/*
 * Copyright 2013 National Bank of Belgium
 * 
 * Licensed under the EUPL, Version 1.1 or – as soon they will be approved 
 * by the European Commission - subsequent versions of the EUPL (the "Licence");
 * You may not use this work except in compliance with the Licence.
 * You may obtain a copy of the Licence at:
 * 
 * http://ec.europa.eu/idabc/eupl
 * 
 * Unless required by applicable law or agreed to in writing, software 
 * distributed under the Licence is distributed on an "AS IS" basis,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the Licence for the specific language governing permissions and 
 * limitations under the Licence.
 */
package jdplus.dfm.base.core;

import java.util.ArrayList;
import java.util.List;
import jdplus.dfm.base.core.var.VarDescriptor;
import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.api.math.Complex;
import jdplus.toolkit.base.api.math.matrices.Matrix;
import jdplus.toolkit.base.core.data.DataBlock;
import jdplus.toolkit.base.core.math.functions.ParamValidation;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;
import jdplus.toolkit.base.core.math.matrices.MatrixException;
import jdplus.toolkit.base.core.math.matrices.decomposition.EigenSystem;
import jdplus.toolkit.base.core.math.matrices.decomposition.IEigenSystem;
import jdplus.toolkit.base.core.ssf.ISsfInitialization;
import jdplus.toolkit.base.core.ssf.multivariate.IMultivariateSsf;
import jdplus.toolkit.base.core.ssf.univariate.ISsfMeasurement;
import lombok.NonNull;

///**
// *
// * @author Jean Palate
// */
public class DfmMapping implements IDfmMapping {

    static final double EPS = 1e-5;

    private final DynamicFactorModel template;
    private final ISsfInitialization.Type initialization;
    private final int nxlags;
    // [0, nml[ loadings
    // [nml, nml+nm[ meas. variance (square roots)
    // [nml+nm, nml+nm+nb*nb*nl[ var parameters 
    // trans. covariance = I
    private final int np;
    private final int nml, nm, nb, nl;
    private final int l0, mv0, v0;
    private final int immax, ifmax;
    private final double cmax;

    private DoubleSeq loadings(DoubleSeq p) {
        return l0 < 0 ? null : p.extract(l0, nml);
    }

    private DoubleSeq vparams(DoubleSeq p) {
        return v0 < 0 ? null : p.extract(v0, nb * nb * nl);
    }

    private DoubleSeq mvars(DoubleSeq p) {
        return mv0 < 0 ? null : p.extract(mv0, nm);
    }

    private DataBlock loadings(DataBlock p) {
        return l0 < 0 ? null : p.extract(l0, nml);
    }

    private DataBlock vparams(DataBlock p) {
        return v0 < 0 ? null : p.extract(v0, nb * nb * nl);
    }

    private DataBlock mvars(DataBlock p) {
        return mv0 < 0 ? null : p.extract(mv0, nm);
    }

    public DfmMapping(DynamicFactorModel model, ISsfInitialization.Type initialization, int nxlags) {
        this(model, initialization, nxlags, false, false);
    }

    public DfmMapping(DynamicFactorModel model, ISsfInitialization.Type initialization, int nxlags, final boolean mfixed, final boolean tfixed) {
        template = model;
        this.initialization=initialization;
        this.nxlags=nxlags;
        nb = template.getNfactors();
        nl = template.getNlags();
        // measurement: all loadings, all var
        // vparams
        // covar
        int p;
        if (mfixed) {
            nml = 0;
            nm = 0;
            l0 = -1;
            mv0 = -1;
            v0 = 0;
            immax = -1;
            ifmax = -1;
            cmax = 0;
            p = nb * nb * nl;
        } else {
            int n = 0, m = 0;
            int im = -1, f = -1;
            double c = 0;
            for (MeasurementDescriptor desc : template.getMeasurements()) {
                for (int i = 0; i < nb; ++i) {
                    double cur = desc.getCoefficient(i);
                    if (!Double.isNaN(cur)) {
                        if (Math.abs(cur) > Math.abs(c)) {
                            c = cur;
                            f = i;
                            im = m;
                        }
                        ++n;
                    }
                }
                ++m;
            }
            l0 = 0;
            immax = im;
            ifmax = f;
            cmax = c;
            nm = template.getMeasurementsCount();
            nml = n - 1;
            mv0 = nml;
            p = nm + nml;
            if (tfixed) {
                v0 = -1;
            } else {
                //         p = tv0 + nb;
                v0 = p;
                p += nb * nb * nl;
            }
        }
        np = p;

    }

    @Override
    public DoubleSeq getDefaultParameters() {
        return map(template);
    }

    @Override
    public IMultivariateSsf map(DoubleSeq p) {
        @NonNull VarDescriptor var = template.getVar();
        DoubleSeq vp = vparams(p);
        if (vp != null) {
            FastMatrix t = FastMatrix.of(var.getCoefficients());
            vp.copyTo(t.getStorage(), 0);
            var=new VarDescriptor(t);
        }
        
        List<MeasurementDescriptor> m = new ArrayList<>();
        
        DoubleSeq l = loadings(p);
        DoubleSeq mv = mvars(p);
        int i0 = 0, j0 = 0;
        if (l != null) {
            int n = 0;
            for (MeasurementDescriptor desc : template.getMeasurements()) {
                MeasurementDescriptor.Builder dbuilder = desc.toBuilder();
                double[] c = desc.getCoefficient().toArray();
                for (int k = 0; k < nb; ++k) {
                    if (!Double.isNaN(c[k])) {
                        if (immax != n || ifmax != k) {
                            c[k] = l.get(i0++);
                        } else {
                            c[k] = cmax;
                        }
                    }
                }
                double x = mv.get(j0++);
                m.add(dbuilder.coefficient(DoubleSeq.of(c))
                        .variance(x * x)
                        .build());
                ++n;
            }
        }
        DynamicFactorModel dfm=new DynamicFactorModel(var, m);
        return dfm.ssfRepresentation(initialization, nxlags);
    }

    @Override
    public DoubleSeq map(DynamicFactorModel m) {
        // copy to p
        DataBlock p = DataBlock.make(np);
        DataBlock l = loadings(p);
        DataBlock mv = mvars(p);
        int i0 = 0, j0 = 0;
        if (l != null) {
            int n = 0;
            for (MeasurementDescriptor desc : m.getMeasurements()) {
                for (int k = 0; k < nb; ++k) {
                    double c=desc.getCoefficient(k);
                    if (!Double.isNaN(c)) {
                        if (n != immax || k != ifmax) {
                            l.set(i0++, c);
                        }
                    }
                }
                mv.set(j0++, Math.sqrt(desc.getVariance()));
                ++n;
            }
        }
        DataBlock vp = vparams(p);
        if (vp != null) {
            Matrix t = m.getVar().getCoefficients();
            vp.copy(t);
        }
        return p;
    }

    @Override
    public boolean checkBoundaries(DoubleSeq inparams) {
        // check the stability of VAR
        try {
            DoubleSeq vp = vparams(inparams);
            if (vp == null) {
                return true;
            }
            // s=(f0,t f1,t f2,t f0,t-1 f1,t-1 f2,t-1 ...f0,t-l+1 f1,t-l+1 f2,t-l+1)
            //    |x00 x10 x20   
            // T =|...
            // T =|1   0   0
            //    |0   1   0
            //    |...
            FastMatrix Q = FastMatrix.square(nb * nl);
            for (int i = 0, i0 = 0; i < nb; ++i) {
                for (int l = 0; l < nl; ++l, i0 += nb) {
                    DataBlock c = Q.column(l * nb + i).range(0, nb);
                    c.copy(vp.extract(i0, nb));
                }
            }
            Q.subDiagonal(-nb).set(1);
            IEigenSystem es = EigenSystem.create(Q, false);
            Complex[] ev = es.getEigenValues();
            for (int i=0; i<ev.length; ++i){
                if (ev[i].absSquare()>=1)
                    return false;
            }
            return true;
        } catch (MatrixException err) {
            return false;
        }
//        return true;
    }
    
    

    @Override
    public double epsilon(DoubleSeq inparams, int idx) {
        return inparams.get(idx) > 0 ? -EPS : EPS;
    }

    @Override
    public int getDim() {
        return np;
    }

    @Override
    public double lbound(int idx) {
        return -Double.MAX_VALUE;
    }

    @Override
    public double ubound(int idx) {
        return Double.MAX_VALUE;
    }

    @Override
    public ParamValidation validate(DataBlock ioparams) {
        return checkBoundaries(ioparams) ? ParamValidation.Valid : ParamValidation.Invalid;
    }

    @Override
    public String getDescription(int idx) {
        return PARAM + idx;
    }
}
