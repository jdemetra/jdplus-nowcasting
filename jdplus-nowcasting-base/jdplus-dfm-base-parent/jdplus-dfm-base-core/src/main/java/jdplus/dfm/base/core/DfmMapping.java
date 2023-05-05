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
import jdplus.toolkit.base.core.math.matrices.SymmetricMatrix;
import jdplus.toolkit.base.core.math.matrices.decomposition.EigenSystem;
import jdplus.toolkit.base.core.math.matrices.decomposition.IEigenSystem;

/**
 *
 * @author palatej
 */
public class DfmMapping implements IDfmMapping {

    private final DynamicFactorModel template;
   // [0, nml[ loadings
    // [nml, nml+nm[ meas. variance (square roots)
    // [nml+nm, nml+nm+nb*nb*nl[ var parameters 
    // [nml+nb*nb*nl+nm, nml+nb*nb*nl+nm+nb*(nb-1)/2[ trans. covariance (cholesky factor), by row 
    private final int np;
    private final int nml, nm, nb, nl;
    private final int l0, mv0, v0, tv0;
    private final int ivmax;
    private final double vmax;
    private final int[] mmax;
    private final double[] fmax;

    private DoubleSeq loadings(DoubleSeq p) {
        return l0 < 0 ? null : p.extract(l0, nml);
    }

    private DoubleSeq vcoefficients(DoubleSeq p) {
        return v0 < 0 ? null : p.extract(v0, nb * nb * nl);
    }

    private DoubleSeq mvars(DoubleSeq p) {
        return mv0 < 0 ? null : p.extract(mv0, nm);
    }

    private DoubleSeq vinnovations(DoubleSeq p) {
        return tv0 < 0 ? null : p.extract(tv0, nb * (nb + 1) / 2);
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

    private DataBlock vinnovations(DataBlock p) {
//        return tv0 < 0 ? null : p.extract(tv0, nb );
        return tv0 < 0 ? null : p.extract(tv0, nb * (nb + 1) / 2);
    }

    private void mvinnovations(FastMatrix v, DoubleSeq vi) {
        int i0 = 0;
        FastMatrix tmp = FastMatrix.square(nb);
        for (int i = 0; i < nb; ++i) {
            DataBlock x = tmp.row(i).range(0, i + 1);
            x.copy(vi.extract(i0, i + 1));
            i0 += i + 1;
        }
        SymmetricMatrix.LLt(tmp, v);
    }

    public DfmMapping(DynamicFactorModel model) {
        this(model, false, false);
    }

    public DfmMapping(DynamicFactorModel model, final boolean fixedMeasurements, final boolean fixedVar) {
        template = model;
        nb = template.getNfactors();
        nl = template.getNlags();
        int p = 0;
        if (fixedMeasurements) {
            nml = 0; // number of loading coefficients
            nm = 0;  // number of measurement var (not fixed)
            l0 = -1; // start of the loading undefined
            mv0 = -1; // start of the measurement variances
            mmax = null; // position of the max measurement coeff 
            fmax = null; // max measurement coeff (fixed)
            ivmax = -1; // position of the max measurement variance
            vmax = 0; // max measurement variance (fixed)
            v0 = 0; // start of the VAR coefficients 
            tv0 = nb * nb * nl; // number of the VAR coefficients
            p = tv0 + nb * (nb + 1) / 2; // number of items in the cholesky factor of the VAR innovations
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
            if (fixedVar) {
                // if the VAR is fixed, we don't fix elements in the measurement (except 1 variance)
                nml = n;
                mmax = null;
                fmax = null;
                mv0 = nml;
                v0 = -1;
                tv0 = -1;
                p = nm + nml;
            } else {
                nml = n - nb;
                mmax = new int[nb];
                for (int i = 0; i < nb; ++i) {
                    mmax[i] = -1;
                }
                n = 0;
                fmax = new double[nb];
                for (MeasurementDescriptor desc : template.getMeasurements()) {
                    for (int j = 0; j < nb; ++j) {
                        double f = desc.getCoefficient(j);
                        if (!Double.isNaN(f) && (mmax[j] < 0 || Math.abs(f) > fmax[j])) {
                            mmax[j] = n;
                            fmax[j] = f;
                        }
                    }
                    ++n;
                }
                for (int i = 0; i < nb; ++i) {
                    if (fmax[i] == 0) {
                        fmax[i] = 1;
                    }
                }
                mv0 = nml;
                //         p = tv0 + nb;
                p += nml + nm;
                v0 = p;
                tv0 = p + nb * nb * nl;
                p = tv0 + nb * (nb + 1) / 2;
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
                        if (mmax == null || n != mmax[k]) {
                            c[k] = lcursor.getAndNext();
                        } else {
                            c[k] = fmax[k];
                        }
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

        DoubleSeq vi = DfmMapping.this.vinnovations(p), vc = DfmMapping.this.vcoefficients(p);
        VarDescriptor var = template.getVar();
        if (vi != null) {
            FastMatrix v = FastMatrix.square(nb);
            mvinnovations(v, vi);
            FastMatrix t = FastMatrix.make(nb, nb * nl);
            vc.copyTo(t.getStorage(), 0);
            var = new VarDescriptor(t, v, var.getInitialization());
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
            DoubleSeqCursor.OnMutable lcursor = l.cursor();
            DoubleSeqCursor.OnMutable vcursor = mv.cursor();
            int n = 0;
            for (MeasurementDescriptor desc : m.getMeasurements()) {
                for (int k = 0; k < nb; ++k) {
                    double c = desc.getCoefficient(k);
                    if (!Double.isNaN(c) && (mmax == null || n != mmax[k])) {
                        lcursor.setAndNext(c);
                    }
                }
                if (n != ivmax) {
                    vcursor.setAndNext(Math.sqrt(desc.getVariance()));
                }
                ++n;
            }
        }
        DataBlock vi = vinnovations(p), vc = vcoefficients(p);
        if (vi != null) {
            FastMatrix v = FastMatrix.of(m.getVar().getInnovationsVariance());
            SymmetricMatrix.lcholesky(v);
            int i0 = 0;
            for (int i = 0; i < nb; ++i) {
                vi.extract(i0, i + 1).copy(v.row(i).range(0, i + 1));
                i0 += i + 1;
            }
            Matrix t = m.getVar().getCoefficients();
            vc.copy(t);
        }
        return p;
    }

    @Override
    public boolean checkBoundaries(DoubleSeq inparams) {
        // check the stability of VAR
        try {
            DoubleSeq vp = vcoefficients(inparams);
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
    public int getDim() {
        return np;
    }

}
