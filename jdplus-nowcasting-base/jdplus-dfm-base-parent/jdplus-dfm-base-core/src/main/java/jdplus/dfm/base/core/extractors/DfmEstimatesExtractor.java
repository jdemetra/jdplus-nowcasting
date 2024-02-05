/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package jdplus.dfm.base.core.extractors;

import jdplus.dfm.base.api.DfmDictionaries;
import jdplus.dfm.base.core.DfmEstimates;
import jdplus.toolkit.base.api.information.InformationExtractor;
import jdplus.toolkit.base.api.information.InformationMapping;
import jdplus.toolkit.base.api.math.matrices.Matrix;
import nbbrd.service.ServiceProvider;

/**
 *
 * @author LEMASSO
 */
@ServiceProvider(InformationExtractor.class)
public class DfmEstimatesExtractor extends InformationMapping<DfmEstimates>{
    
    public DfmEstimatesExtractor() {
        set(DfmDictionaries.LIKELIHOOD_LL, Double.class, source -> source.getLl());
        set(DfmDictionaries.HESSIAN, Matrix.class, source -> source.getHessian());
        set(DfmDictionaries.GRADIENT, double[].class, source -> source.getGradient());
        set(DfmDictionaries.HAS_CONVERGED, Boolean.class, source -> source.isHasConverged());
    }
    
    @Override
    public Class<DfmEstimates> getSourceClass() {
        return DfmEstimates.class;
    }
}
