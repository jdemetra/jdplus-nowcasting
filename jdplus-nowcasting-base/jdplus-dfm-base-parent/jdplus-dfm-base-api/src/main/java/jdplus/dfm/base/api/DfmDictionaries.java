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
            INPUT_TRANSFORMED = "input_transformed", MEASUREMENT_COEFFICIENTS = "measurement_coefficients",
            MEASUREMENT_ERRORS_VARIANCE = "measurement_errors_variance", VAR_COEFFICIENTS = "var_coefficients", 
            VAR_ERRORS_VARIANCE = "var_errors_variance", FACTORS = "factors", FACTORS_STDERR = "factors_stderr", 
            FORECASTS_TRANSFORMED="forecasts_transformed", FORECASTS_TRANSFORMED_STDERR = "forecasts_transformed_stderr", 
            FORECASTS="forecasts", FORECASTS_STDERR = "forecasts_stderr", RESIDUALS = "residuals", 
            RESIDUALS_STANDARDIZED = "residuals_standardized", RESIDUALS_STANDARDIZED_SMOOTHED = "residuals_standardized_smoothed", 
            LIKELIHOOD_LL = "likelihood_ll", INITIALIZATION_TYPE = "initialization_type", FACTORS_TYPE = "factors_type",
            HESSIAN = "hessian", GRADIENT = "gradient", HAS_CONVERGED = "has_converged";
            
    public final Dictionary DFMDICTIONARY = AtomicDictionary.builder()
            .name("dfm")
            .item(AtomicDictionary.Item.builder().name(SAMPLE_MEAN).description("Sample mean of the original series").outputClass(double[].class).build())
            .item(AtomicDictionary.Item.builder().name(SAMPLE_STDDEV).description("Sample standard deviation of the original series").outputClass(double[].class).build())
            .item(AtomicDictionary.Item.builder().name(INPUT).description("input data").outputClass(Matrix.class).build()) 
            .item(AtomicDictionary.Item.builder().name(INPUT_TRANSFORMED).description("Standardized input data, corrected for mean and stdev").outputClass(Matrix.class).build()) 
            .item(AtomicDictionary.Item.builder().name(MEASUREMENT_COEFFICIENTS).description("Estimated coefficients related to the normalized factors in the measurement equation").outputClass(Matrix.class).build())
            .item(AtomicDictionary.Item.builder().name(MEASUREMENT_ERRORS_VARIANCE).description("Variance of the error terms in the measurement equation").outputClass(double[].class).build())
            .item(AtomicDictionary.Item.builder().name(VAR_COEFFICIENTS).description("Estimated coefficients related to the Vector Autoregressive Process").outputClass(Matrix.class).build())
            .item(AtomicDictionary.Item.builder().name(VAR_ERRORS_VARIANCE).description("Matrice variance-covariance of the error terms related to the Vector Autoregressive Process").outputClass(Matrix.class).build())
            .item(AtomicDictionary.Item.builder().name(FACTORS).description("Estimated factors").outputClass(Matrix.class).build())
            .item(AtomicDictionary.Item.builder().name(FACTORS_STDERR).description("Standard error of the estimated factors").outputClass(Matrix.class).build())
            .item(AtomicDictionary.Item.builder().name(FORECASTS_TRANSFORMED).description("Forecasts of the transformed input over a one-year horizon").outputClass(TsData.class).build())
            .item(AtomicDictionary.Item.builder().name(FORECASTS_TRANSFORMED_STDERR).description("Standard error of the forecasts of the transformed input").outputClass(TsData.class).build())
            .item(AtomicDictionary.Item.builder().name(FORECASTS).description("Forecasts over a one-year horizon").outputClass(TsData.class).build())
            .item(AtomicDictionary.Item.builder().name(FORECASTS_STDERR).description("Standard error of the forecasts").outputClass(TsData.class).build())
            .item(AtomicDictionary.Item.builder().name(RESIDUALS).description("One-step forecast errors").outputClass(Matrix.class).build())
            .item(AtomicDictionary.Item.builder().name(RESIDUALS_STANDARDIZED).description("Standardized one-step forecast errors").outputClass(Matrix.class).build())
            .item(AtomicDictionary.Item.builder().name(LIKELIHOOD_LL).description("Log-likelihood").outputClass(Double.class).build())
            .item(AtomicDictionary.Item.builder().name(INITIALIZATION_TYPE).description("Type of initialization").outputClass(String.class).build())
            .item(AtomicDictionary.Item.builder().name(FACTORS_TYPE).description("Type of each factors, represented by the number of lags implied by the measurement").outputClass(int[].class).build())
            .item(AtomicDictionary.Item.builder().name(HESSIAN).description("HESSIAN matrix related to ML estimation of the parameters").outputClass(Matrix.class).build())
            .item(AtomicDictionary.Item.builder().name(GRADIENT).description("GRANDIENT vector related to ML estimation of the parameters").outputClass(double[].class).build())
            .item(AtomicDictionary.Item.builder().name(HAS_CONVERGED).description("Indication whether ML estimation of the parameters has converged").outputClass(Boolean.class).build())
            .build();
}
