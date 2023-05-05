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

import jdplus.toolkit.base.api.data.DoubleSeq;
import jdplus.toolkit.base.api.math.matrices.Matrix;
import jdplus.toolkit.base.core.data.DataBlock;
import jdplus.toolkit.base.core.math.linearsystem.LinearSystemSolver;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;
import jdplus.toolkit.base.core.math.matrices.SymmetricMatrix;
import jdplus.toolkit.base.core.ssf.ISsfInitialization;

/**
 *
 * @author Jean Palate
 */
public class Initialization implements ISsfInitialization {

    public static Initialization zero(Matrix V, int nxlags) {
        int nf=V.getRowsCount();
        int dim = nf * nxlags;
        FastMatrix V0 = FastMatrix.square(dim);
        for (int i = 0; i < nf; ++i) {
            DataBlock cv = V0.column(i * nxlags).extract(0, nf, nxlags);
            cv.set(V.column(i));
        }
        return new Initialization(null, V0, dim);
    }

    public static Initialization unconditional(Dynamics dynamics) {
        FastMatrix v0 = of(dynamics);
        return new Initialization(null, v0, v0.getRowsCount());
    }

    public static Initialization user(DoubleSeq a0, Matrix V0) {
        return new Initialization(a0, FastMatrix.of(V0), V0.getRowsCount());
    }

    public static FastMatrix of(Dynamics dynamics) {
        int nl = dynamics.nl(), nf = dynamics.nf();
        // We have to solve the steady state equation:
        // V = T V T' + Q
        // We consider first the [nl*nf, nl*nf] sub-system
        FastMatrix v = dynamics.getV();
        FastMatrix t = dynamics.getT();

        int n = nf * nl;
        FastMatrix cov = FastMatrix.square(n);
        int np = (n * (n + 1)) / 2;
        FastMatrix M = FastMatrix.square(np);
        double[] b = new double[np];
        // fill the matrix
        for (int c = 0, i = 0; c < n; ++c) {
            for (int r = c; r < n; ++r, ++i) {
                M.set(i, i, 1);
                if (r % nl == 0 && c % nl == 0) {
                    b[i] = v.get(r / nl, c / nl);
                }
                for (int k = 0; k < n; ++k) {
                    for (int l = 0; l < n; ++l) {
                        double zr = 0, zc = 0;
                        if (r % nl == 0) {
                            zr = t.get(r / nl, l);
                        } else if (r == l + 1) {
                            zr = 1;
                        }
                        if (c % nl == 0) {
                            zc = t.get(c / nl, k);
                        } else if (c == k + 1) {
                            zc = 1;
                        }
                        double z = zr * zc;
                        if (z != 0) {
                            int p = l <= k ? pos(k, l, n) : pos(l, k, n);
                            M.add(i, p, -z);
                        }
                    }
                }
            }
        }

        LinearSystemSolver.robustSolver().solve(M, DataBlock.of(b));

        for (int i = 0, j = 0; i < n; i++) {
            cov.column(i).drop(i, 0).copyFrom(b, j);
            j += n - i;
        }
        SymmetricMatrix.fromLower(cov);
        int nlx = dynamics.nxlags;
        if (nl == nlx) {
            return cov;
        }
        int dim = nlx * nf;
        FastMatrix fullCov = FastMatrix.square(dim);

        for (int r = 0; r < nf; ++r) {
            for (int c = 0; c < nf; ++c) {
                fullCov.extract(r * nlx, nl, c * nlx, nl).copy(cov.extract(r * nl, nl, c * nl, nl));
            }
        }
        for (int i = nl; i < nlx; ++i) {
            dynamics.TVT(0, fullCov);
            dynamics.addV(0, fullCov);
        }
        return fullCov;
    }

    private Initialization(DoubleSeq a0, FastMatrix V0, int dim) {
        this.a0 = a0;
        this.V0 = V0;
        this.dim = dim;
    }

    private final DoubleSeq a0;
    private final int dim;
    private final FastMatrix V0;

    @Override
    public int getStateDim() {
        return dim;
    }

    @Override
    public boolean isDiffuse() {
        return false;
    }

    @Override
    public int getDiffuseDim() {
        return 0;
    }

    @Override
    public void diffuseConstraints(FastMatrix b) {
    }

    @Override
    public void a0(DataBlock a0) {
        if (this.a0 != null) {
            a0.copy(this.a0);
        }
    }

    @Override
    public void Pf0(FastMatrix pf0) {
        if (V0 != null) {
            pf0.copy(V0);
        }
    }

    private static int pos(int r, int c, int n) {
        return r + c * (2 * n - c - 1) / 2;
    }
}
