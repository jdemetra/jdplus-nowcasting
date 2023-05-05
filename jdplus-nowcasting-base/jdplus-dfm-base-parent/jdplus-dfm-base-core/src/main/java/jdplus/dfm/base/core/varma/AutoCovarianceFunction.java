/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */
package jdplus.dfm.base.core.varma;

import java.util.ArrayList;
import jdplus.toolkit.base.core.math.matrices.FastMatrix;
import jdplus.toolkit.base.core.math.matrices.GeneralMatrix;

/**
 *
 * @author palatej
 */
public class AutoCovarianceFunction {

    private final ArrayList<FastMatrix> covs = new ArrayList<>();
    private final ArrayList<FastMatrix> g = new ArrayList<>();
    private final VarmaModel varma;

    public AutoCovarianceFunction(VarmaModel varma) {
        this.varma = varma;
        g.add(varma.sig());
        computeInitialCov();

    }
    
    public FastMatrix cov(int lag) {
        
        if (lag >= covs.size()) {
            compute(lag);
        }
        return covs.get(lag);
    }
    
    private void computeInitialCov(){
        
    }

    private void compute(int lag) {
        computeg(lag);
        for (int i=covs.size(); i<=lag; ++i){
            calcCov(i);
        }
    }

    private void computeg(int lag) {
        int n=varma.getDim(), p=varma.getP(), q=varma.getQ();
        if (g.size() > lag) {
            return;
        }
        for (int i = g.size(); i <= lag; ++i) {
            FastMatrix g = FastMatrix.square(n);
            if (i <= q) {
                g = GeneralMatrix.AB(varma.th(i),varma.sig());
            }
            for (int j = 1; j <= Math.min(i, p); ++j) {
                FastMatrix tmp = GeneralMatrix.AB(varma.phi(j), this.g.get(i-j));
                g.sub(tmp);
            }
            this.g.add(g);
        }
    }

    private void calcCov(int k) {
        int n=varma.getDim(), p=varma.getP(), q=varma.getQ();
        FastMatrix cov=g.get(k).deepClone();
        for (int i=1; i<=q; ++i){
            cov.add(GeneralMatrix.AB(g.get(k-i), varma.th(k-i)));
        }
    }

}
