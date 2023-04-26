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

import jdplus.toolkit.base.core.data.DataBlock;
import jdplus.toolkit.base.core.data.DataWindow;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;
import jdplus.toolkit.base.core.ssf.ISsfDynamics;

/**
 *
 * @author Jean Palate
 */
public class Dynamics implements ISsfDynamics {

    public static Dynamics of(FastMatrix T, FastMatrix V) {
        return new Dynamics(T, V);
    }

    private Dynamics(FastMatrix T, FastMatrix V, int nxlags) {
        this.T = T;
        this.V = V;
        int nf = V.getRowsCount();
        ttmp = new double[nf];
        xtmp = new double[nf * nxlags];
        this.nxlags=nxlags;
    }

    public static Dynamics of(FastMatrix T, FastMatrix V, int nxlags) {
        return new Dynamics(T, V, nxlags);
    }

    private Dynamics(FastMatrix T, FastMatrix V) {
        this.T = T;
        this.V = V;
        ttmp = new double[V.getColumnsCount()];
        xtmp = new double[T.getColumnsCount()];
        nxlags=xtmp.length/ttmp.length;
    }

    final int nxlags;
    final FastMatrix T, V;
    final double[] ttmp, xtmp;
    
    public FastMatrix getT(){
        return T;
    }
    
    public FastMatrix getV(){
        return V;
    }

    int nl() {
        return T.getColumnsCount() / T.getRowsCount();
    }

    int nf() {
        return T.getRowsCount();
    }

    @Override
    public int getInnovationsDim() {
        return V.getRowsCount();
    }

    @Override
    public void V(int pos, FastMatrix qm) {
        qm.copy(V);
    }

    @Override
    public void S(int pos, FastMatrix cm) {
        int nf = nf();
        for (int i = 0, r = 0; i < nf; ++i, r += nxlags) {
            cm.set(r, i, 1);
        }
    }

    @Override
    public boolean hasInnovations(int pos) {
        return true;
    }

    @Override
    public boolean areInnovationsTimeInvariant() {
        return true;
    }

    @Override
    public void T(int pos, FastMatrix tr) {
        int nl = nl(), nf = nf();
        for (int i = 0, r = 0; i < nf; ++i, r += nxlags) {
            for (int j = 0, c = 0; j < nf; ++j, c += nxlags) {
                FastMatrix B = tr.extract(r, r + nxlags, c, c + nxlags);
                if (i == j) {
                    B.subDiagonal(-1).set(1);
                }
                B.row(0).range(0, nl).copy(T.row(i).range(j * nl, (j + 1) * nl));
            }
        }
    }

    @Override
    public void TX(int pos, DataBlock x) {
        int nl = nl(), nf = nf();
        // compute first the next item
        for (int i = 0; i < nf; ++i) {
            DataWindow p = T.row(i).left();
            DataWindow xb = x.left();
            double r = p.next(nl).dot(xb.next(nl));
            for (int j = 1; j < nf; ++j) {
                r += p.next(nl).dot(xb.move(nxlags));
            }
            ttmp[i] = r;
        }
        x.fshiftAndZero();
        x.extract(0, nf, nxlags).copyFrom(ttmp, 0);
    }

    @Override
    public void addSU(int pos, DataBlock x, DataBlock u) {
        x.extract(0, nf(), nxlags).add(u);
    }

    @Override
    public void addV(int pos, FastMatrix p) {
        int nf = nf();
        for (int i = 0; i < nf; ++i) {
            DataBlock cv = p.column(i * nxlags).extract(0, nf, nxlags);
            cv.add(V.column(i));
        }
    }

    @Override
    public void XT(int pos, DataBlock x) {
        int nl = nl(), nf = nf();
        for (int i = 0, k = 0, l = 0; i < nf; ++i) {
            for (int j = 0; j < nl; ++j, ++k) {
                double r = ((k + 1) % nxlags != 0) ? x.get(k + 1) : 0;
                r += T.column(l++).dot(x.extract(0, nf, nxlags));
                xtmp[k] = r;
            }
            for (int j = nl; j < nxlags - 1; ++j, ++k) {
                xtmp[k] = x.get(k + 1);
            }
            if (nxlags > nl) {
                xtmp[k++] = 0;
            }
        }
        x.copyFrom(xtmp, 0);
    }

    @Override
    public void XS(int pos, DataBlock x, DataBlock xs) {
        xs.copy(x.extract(0, nf(), nxlags));
    }

    @Override
    public boolean isTimeInvariant() {
        return true;
    }
}
