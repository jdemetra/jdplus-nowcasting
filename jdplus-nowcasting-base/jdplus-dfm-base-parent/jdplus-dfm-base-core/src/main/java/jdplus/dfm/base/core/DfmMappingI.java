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

import java.util.ArrayList;
import java.util.List;
import jdplus.dfm.base.core.var.VarDescriptor;
import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.api.data.DoubleSeqCursor;
import jdplus.toolkit.base.api.math.Complex;
import jdplus.toolkit.base.api.math.matrices.Matrix;
import jdplus.toolkit.base.core.data.DataBlock;
import jdplus.toolkit.base.core.math.functions.ParamValidation;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;
import jdplus.toolkit.base.core.math.matrices.MatrixException;
import jdplus.toolkit.base.core.math.matrices.decomposition.EigenSystem;
import jdplus.toolkit.base.core.math.matrices.decomposition.IEigenSystem;

///**
// * Mapping for models with independent shocks in the VAR
// * @author Jean Palate
// */
public class DfmMappingI implements IDfmMapping {

    private final DynamicFactorModel template;
    // [0, nml[ loadings
    // [nml, nml+nm[ meas. variance (square roots)
    // [nml+nm, nml+nm+nb*nb*nl[ var parameters 
    // trans. covariance = I
    private final int np;
    private final int nml, nm, nb, nl;
    private final int l0, mv0, v0;
    private final int ivmax;
    private final double vmax;

    private DoubleSeq loadings(DoubleSeq p) {
        return l0 < 0 ? null : p.extract(l0, nml);
    }

    private DoubleSeq vcoefficients(DoubleSeq p) {
        return v0 < 0 ? null : p.extract(v0, nb * nb * nl);
    }

    private DoubleSeq mvars(DoubleSeq p) {
        return mv0 < 0 ? null : p.extract(mv0, nm);
    }

    private DataBlock loadings(DataBlock p) {
        return l0 < 0 ? null : p.extract(l0, nml);
    }

    private DataBlock vcoefficients(DataBlock p) {
        return v0 < 0 ? null : p.extract(v0, nb * nb * nl);
    }

    private DataBlock mvars(DataBlock p) {
        return mv0 < 0 ? null : p.extract(mv0, nm);
    }

    public DfmMappingI(DynamicFactorModel model) {
        this(model, false, false);
    }

    public DfmMappingI(DynamicFactorModel model, final boolean mfixed, final boolean tfixed) {
        template = model.normalize();
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
            ivmax = -1;
            vmax = 0;
            p = nb * nb * nl;
        } else {
            int n = 0, m = 0;
            int iv = -1;
            double v = 0;
            for (MeasurementDescriptor desc : template.getMeasurements()) {
                if (desc.getVariance() > v) {
                    v = desc.getVariance();
                    iv = m;
                }
                for (int i = 0; i < nb; ++i) {
                    double cur = desc.getCoefficient(i);
                    if (!Double.isNaN(cur)) {
                        ++n;
                    }
                }
                ++m;
            }
            l0 = 0;
            ivmax = iv;
            vmax = v;
            nm = template.getMeasurementsCount() - 1;
            nml = n;
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
    public DynamicFactorModel map(DoubleSeq p) {
        VarDescriptor var = template.getVar();
        DoubleSeq vp = vcoefficients(p);
        if (vp != null) {
            FastMatrix t = FastMatrix.of(var.getCoefficients());
            vp.copyTo(t.getStorage(), 0);
            var = new VarDescriptor(t);
        }
        List<MeasurementDescriptor> m;
        DoubleSeq l = loadings(p);
        DoubleSeq mv = mvars(p);
        if (l != null) {
            m = new ArrayList<>();
            DoubleSeqCursor lcursor = l.cursor();
            DoubleSeqCursor vcursor = mv.cursor();
            int n = 0;
            for (MeasurementDescriptor desc : template.getMeasurements()) {
                MeasurementDescriptor.Builder dbuilder = desc.toBuilder();
                double[] c = desc.getCoefficient().toArray();
                for (int k = 0; k < nb; ++k) {
                    if (!Double.isNaN(c[k])) {
                        c[k] = lcursor.getAndNext();
                    }
                }
                if (n != ivmax) {
                    double x = vcursor.getAndNext();
                    m.add(dbuilder.coefficient(DoubleSeq.of(c))
                            .variance(x * x)
                            .build());
                } else {
                    m.add(dbuilder.coefficient(DoubleSeq.of(c))
                            .variance(vmax)
                            .build());
                }
                ++n;
            }
        } else {
            m = template.getMeasurements();
        }
        return new DynamicFactorModel(var, m);
    }

    @Override
    public DoubleSeq map(DynamicFactorModel m) {
        // copy to p
        DataBlock p = DataBlock.make(np);
        DataBlock l = loadings(p);
        DataBlock mv = mvars(p);
        if (l != null) {
            int n = 0;
            DoubleSeqCursor.OnMutable lcursor = l.cursor();
            DoubleSeqCursor.OnMutable vcursor = mv.cursor();
            for (MeasurementDescriptor desc : m.getMeasurements()) {
                for (int k = 0; k < nb; ++k) {
                    double c = desc.getCoefficient(k);
                    if (!Double.isNaN(c)) {
                        lcursor.setAndNext(c);
                    }
                }
                if (n != ivmax) {
                    vcursor.setAndNext(Math.sqrt(desc.getVariance()));
                }
                ++n;
            }
        }
        DataBlock vc = vcoefficients(p);
        if (vc != null) {
            Matrix t = m.getVar().getCoefficients();
            vc.copy(t);
        }
        return p;
    }

    @Override
    public boolean checkBoundaries(DoubleSeq inparams) {
        // check the stability of VAR
        try {
            DoubleSeq vp = DfmMappingI.this.vcoefficients(inparams);
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
            for (int i = 0, i0 = 0; i < Q.getColumnsCount(); ++i, i0 += nb) {
                DataBlock c = Q.column(i).range(0, nb);
                c.copy(vp.extract(i0, nb));
            }

            Q.subDiagonal(-nb).set(1);
            IEigenSystem es = EigenSystem.create(Q, false);
            Complex[] ev = es.getEigenValues();
            for (int i = 0; i < ev.length; ++i) {
                if (ev[i].absSquare() >= 1) {
                    return false;
                }
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
