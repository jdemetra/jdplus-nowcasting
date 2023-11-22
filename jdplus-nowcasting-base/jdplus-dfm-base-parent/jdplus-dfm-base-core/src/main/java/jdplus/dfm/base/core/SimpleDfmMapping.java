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
import jdplus.toolkit.base.api.math.matrices.Matrix;
import jdplus.toolkit.base.core.data.DataBlock;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;

/**
 * Mapping defined for a simplified model (Innovation variance = I, nlags = 1)
 * We just fix the max variance of the measurements
 *
 * @author Jean Palate
 */
public class SimpleDfmMapping implements IDfmMapping {
    
    private final DynamicFactorModel template;
    // [0, nml[ loadings
    // [nml, nml+nm[ meas. variance (square roots). The max is fixed, which means that nm = number of measurements - 1
    // [nml+nm, nml+nm+nb*nb*nl[ var parameters 
    private final int ivmax;
//    private final double vmax;
    private final int np;
    private final int nml, nm, nb, nl;
    private final int v0;
    
    private DoubleSeq loadings(DoubleSeq p) {
        return p.extract(0, nml);
    }
    
    private DoubleSeq vcoefficients(DoubleSeq p) {
        return v0 < 0 ? null : p.extract(v0, nb * nb * nl);
    }
    
    private DoubleSeq mvars(DoubleSeq p) {
        return p.extract(nml, nm);
    }
    
    private DataBlock loadings(DataBlock p) {
        return p.extract(0, nml);
    }
    
    private DataBlock vcoefficients(DataBlock p) {
        return v0 < 0 ? null : p.extract(v0, nb * nb * nl);
    }
    
    private DataBlock mvars(DataBlock p) {
        return p.extract(nml, nm);
    }
    
    public SimpleDfmMapping(DynamicFactorModel model) {
        DynamicFactorModel nmodel = model.normalize();
        nb = model.getNfactors();
        nl = 1;
        FastMatrix C = FastMatrix.of(nmodel.getVar().getCoefficients());
        for (int c = nb; c < C.getColumnsCount(); ++c) {
            C.column(c).set(0);
        }
        template = new DynamicFactorModel(new VarDescriptor(C, model.getVar().getInitialization()), nmodel.getMeasurements());
        nm = nmodel.getMeasurementsCount() - 1;
        // measurement: all loadings, all var
        // vparams
        // covar
        int p = 0;
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
        ivmax = iv;
        nml = n; // number of loadings
        p += nml + nm;
        v0 = p;
        p += nl * nb * nb;
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
                            .variance(desc.getVariance())
                            .build());
                    
                }
                ++n;
            }
        } else {
            m = template.getMeasurements();
        }
        
        DoubleSeq vc = vcoefficients(p);
        FastMatrix t = FastMatrix.make(template.getNfactors(), template.getNfactors()*template.getNlags());
        vc.copyTo(t.getStorage(), 0);
        return new DynamicFactorModel(new VarDescriptor(t, template.getVar().getInitialization()), m);
    }
    
    @Override
    public DoubleSeq map(DynamicFactorModel m) {
        // copy to p
        DataBlock p = DataBlock.make(np);
        DataBlock l = loadings(p);
        DataBlock mv = mvars(p);
        DoubleSeqCursor.OnMutable lcursor = l.cursor();
        DoubleSeqCursor.OnMutable vcursor = mv.cursor();
        int n = 0;
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
        DataBlock vc = vcoefficients(p);
        Matrix t = m.getVar().getCoefficients();
        // just keep the first lag
        for (int k = 0, c=0; k < nb; ++k, c+=nb) {
            vc.extract(c, nb).copy(t.column(k));
        }
        return p;
    }
    
    @Override
    public boolean checkBoundaries(DoubleSeq inparams) {
       DoubleSeq vp = vcoefficients(inparams);
        if (vp != null) {
            int i0 = 0;
            for (int i = 0; i < nb; ++i) {
                for (int j = 0; j < nb; ++j) {
                    if (Math.abs(vp.get(i0++)) > .99) {
                        return false;
                    }
                }
            }
        }
        return true;
    }
    
    @Override
    public int getDim() {
        return np;
    }
    
    public DynamicFactorModel validate(DynamicFactorModel model) {
        Matrix m = model.getVar().getCoefficients();
        FastMatrix v = FastMatrix.make(m.getRowsCount(), m.getColumnsCount());
        DoubleSeq md = m.diagonal();
        DataBlock vd = v.diagonal();
        for (int i = 0; i < nb; ++i) {
            double r = md.get(i);
            if (Math.abs(r) > 1) {
                r = Math.signum(r) * Math.min(.99, 1 / Math.abs(r));
            }
            vd.set(i, r);
        }
        return new DynamicFactorModel(model.getVar().withCoefficients(v), model.getMeasurements());
    }

}
