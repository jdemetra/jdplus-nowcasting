/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package jdplus.dfm.base.core.extractors;

import jdplus.dfm.base.api.DfmDictionaries;
import jdplus.dfm.base.core.DfmResults;
import jdplus.dfm.base.core.DynamicFactorModel;
import jdplus.toolkit.base.api.information.InformationExtractor;
import jdplus.toolkit.base.api.information.InformationMapping;
import jdplus.toolkit.base.api.math.matrices.Matrix;
import nbbrd.service.ServiceProvider;

/**
 *
 * @author LEMASSO
 */
@ServiceProvider(InformationExtractor.class)
public class DfmResultsExtractor extends InformationMapping<DfmResults> {  
   public final int NFCAST = -1;
   
    public DfmResultsExtractor() {
        set(DfmDictionaries.SAMPLE_MEAN, double[].class, source -> source.getSampleMean().toArray());
        set(DfmDictionaries.SAMPLE_STDDEV, double[].class, source -> source.getSampleStDev().toArray());
        set(DfmDictionaries.INPUT, Matrix.class, source -> source.getInputData());
        set(DfmDictionaries.INPUT_TRANSFORMED, Matrix.class, source -> source.getTransformedInputData());
        set(DfmDictionaries.FACTORS, Matrix.class, source -> source.getFactors());
        set(DfmDictionaries.FACTORS_STDERR, Matrix.class, source -> source.getFactorsStdErr());
        set(DfmDictionaries.RESIDUALS, Matrix.class, source -> source.getResiduals());
        set(DfmDictionaries.RESIDUALS_STANDARDIZED, Matrix.class, source -> source.getResidualsStandardized());
        set(DfmDictionaries.LIKELIHOOD_LL, Double.class, source -> source.getLogLikelihood());
        
        delegate(null, DynamicFactorModel.class, source -> source.getDfm());
        
        setArray(DfmDictionaries.FORECASTS_TRANSFORMED, NFCAST, Matrix.class, (source, i) -> source.forecastsT(i));
        setArray(DfmDictionaries.FORECASTS, NFCAST, Matrix.class, (source, i) -> source.forecasts(i));
        setArray(DfmDictionaries.FORECASTS_TRANSFORMED_STDERR, NFCAST, Matrix.class, (source, i) -> source.forecastsTStDev(i));
        setArray(DfmDictionaries.FORECASTS_STDERR, NFCAST, Matrix.class, (source, i) -> source.forecastsStDev(i));      
    }
    
    @Override
    public Class<DfmResults> getSourceClass() {
        return DfmResults.class;
    }
}
    
    

