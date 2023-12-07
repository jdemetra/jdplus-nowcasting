/*
 * Click nbfs://nbhost/SystemFileSystem/Templates/Licenses/license-default.txt to change this license
 * Click nbfs://nbhost/SystemFileSystem/Templates/Classes/Class.java to edit this template
 */
package jdplus.dfm.base.api;

import jdplus.toolkit.base.api.dictionaries.AtomicDictionary;
import jdplus.toolkit.base.api.dictionaries.Dictionary;
import jdplus.toolkit.base.api.timeseries.TsData;
import jdplus.toolkit.base.api.math.matrices.Matrix;
import nbbrd.design.Development;

/**
 *
 * @author LEMASSO
 */
@lombok.experimental.UtilityClass
@Development(status = Development.Status.Beta)
public class DfmDictionaries {
    
    public static final String SAMPLE_MEAN = "sample_mean", SAMPLE_STDDEV = "sample_stddev", INPUT = "input",
            INPUT_TRANSFORMED = "input_transformed", PARAMETERS_FACTORS = "parameters_factors",
            PARAMETERS_FACTORS_VARIANCE = "parameters_factors_variance", PARAMETERS_VAR = "parameters_var", 
            PARAMETERS_VAR_VARIANCE = "parameters_var_variance", FACTORS = "factors", FACTORS_STDERR = "factors_stderr", 
            FORECASTS_TRANSFORMED="forecasts_transformed", FORECASTS_TRANSFORMED_STDERR = "forecasts_transformed_stderr", 
            FORECASTS="forecasts", FORECASTS_STDERR = "forecasts_stderr", RESIDUALS = "residuals", 
            RESIDUALS_STANDARDIZED = "residuals_standardized", RESIDUALS_STANDARDIZED_SMOOTHED = "residuals_standardized_smoothed", 
            LIKELIHOOD_LL = "likelihood_ll";
            
    public final Dictionary DFMDICTIONARY = AtomicDictionary.builder()
            .name("dfm")
            .item(AtomicDictionary.Item.builder().name(SAMPLE_MEAN).description("Sample mean of the original series").outputClass(double[].class).build())
            .item(AtomicDictionary.Item.builder().name(SAMPLE_STDDEV).description("Sample standard deviation of the original series").outputClass(double[].class).build())
            .item(AtomicDictionary.Item.builder().name(INPUT).description("input data").outputClass(Matrix.class).build()) 
            .item(AtomicDictionary.Item.builder().name(INPUT_TRANSFORMED).description("Standardized input data, corrected for mean and stdev").outputClass(Matrix.class).build()) 
            .item(AtomicDictionary.Item.builder().name(PARAMETERS_FACTORS).description("Estimated parameters related to the normalized factors").outputClass(Matrix.class).build())
            .item(AtomicDictionary.Item.builder().name(PARAMETERS_FACTORS_VARIANCE).description("Idiosyncratic variance of the parameters related to the normalized factors").outputClass(double[].class).build())
            .item(AtomicDictionary.Item.builder().name(PARAMETERS_VAR).description("Estimated parameters related to the vector autoregressive process").outputClass(Matrix.class).build())
            .item(AtomicDictionary.Item.builder().name(PARAMETERS_VAR_VARIANCE).description("Matrice variance-covariance of the parameters related to the vector autoregressive process").outputClass(Matrix.class).build())
            .item(AtomicDictionary.Item.builder().name(FACTORS).description("Estimated factors").outputClass(Matrix.class).build())
            .item(AtomicDictionary.Item.builder().name(FACTORS_STDERR).description("Standard error of the estimated factors").outputClass(Matrix.class).build())
            .item(AtomicDictionary.Item.builder().name(FORECASTS_TRANSFORMED).description("Forecasts of the transformed input over a one-year horizon").outputClass(TsData.class).build())
            .item(AtomicDictionary.Item.builder().name(FORECASTS_TRANSFORMED_STDERR).description("Standard error of the forecasts of the transformed input").outputClass(TsData.class).build())
            .item(AtomicDictionary.Item.builder().name(FORECASTS).description("Forecasts over a one-year horizon").outputClass(TsData.class).build())
            .item(AtomicDictionary.Item.builder().name(FORECASTS_STDERR).description("Standard error of the forecasts").outputClass(TsData.class).build())
            .item(AtomicDictionary.Item.builder().name(RESIDUALS).description("One-step forecast errors").outputClass(Matrix.class).build())
            .item(AtomicDictionary.Item.builder().name(RESIDUALS_STANDARDIZED).description("Standardized one-step forecast errors").outputClass(Matrix.class).build())
            .item(AtomicDictionary.Item.builder().name(LIKELIHOOD_LL).description("Log-likelihood").outputClass(Double.class).build())
            .build();
}
