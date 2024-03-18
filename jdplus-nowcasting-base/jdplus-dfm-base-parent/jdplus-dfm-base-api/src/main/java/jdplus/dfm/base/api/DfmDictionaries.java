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
    
    public static final String SERIES_INDEX = "series_index", SERIES_PERIOD = "series_period", 
            SERIES_EXPECTED_VALUE_T = "series_expected_value_T", SERIES_EXPECTED_VALUE = "series_expected_value",
            SERIES_OBSERVED_VALUE_T = "series_observed_value_T", SERIES_OBSERVED_VALUE = "series_observed_value",
            SERIES_NEWS_T = "series_news_T", SERIES_NEWS = "series_news",
            SERIES_WEIGHTS_T = "series_weights_T", SERIES_WEIGHTS = "series_weights",
            SERIES_IMPACTS_T = "series_impacts_T", SERIES_IMPACTS = "series_impacts", FORECASTS_PERIODS = "forecasts_periods",
            OLD_FORECASTS_T = "old_forecasts_T", OLD_FORECASTS = "old_forecasts",
            REVISED_FORECASTS_T = "revised_forecasts_T", REVISED_FORECASTS = "revised_forecasts",
            NEW_FORECASTS_T = "new_forecasts_T", NEW_FORECASTS = "new_forecasts";
    
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
            .item(AtomicDictionary.Item.builder().name(HESSIAN).description("Hessian matrix related to ML estimation of the parameters").outputClass(Matrix.class).build())
            .item(AtomicDictionary.Item.builder().name(GRADIENT).description("Gradient vector related to ML estimation of the parameters").outputClass(double[].class).build())
            .item(AtomicDictionary.Item.builder().name(HAS_CONVERGED).description("Indication whether ML estimation of the parameters has converged").outputClass(Boolean.class).build())
            .build();
    
    public final Dictionary DFMNEWSDICTIONARY = AtomicDictionary.builder()
        .name("dfmNews")
        .item(AtomicDictionary.Item.builder().name(SERIES_INDEX).description("Index number of series where new data was added").outputClass(int[].class).build())
        .item(AtomicDictionary.Item.builder().name(SERIES_PERIOD).description("Period where new data was added").outputClass(String[].class).build())
        .item(AtomicDictionary.Item.builder().name(SERIES_EXPECTED_VALUE_T).description("Previous forecast of the standardized series for the period where new data was added").outputClass(double[].class).build())
        .item(AtomicDictionary.Item.builder().name(SERIES_EXPECTED_VALUE).description("Previous forecast of the original series for the period where new data was added").outputClass(double[].class).build())
        .item(AtomicDictionary.Item.builder().name(SERIES_OBSERVED_VALUE_T).description("Observed value of the standardized series for the period where new data was added").outputClass(double[].class).build())
        .item(AtomicDictionary.Item.builder().name(SERIES_OBSERVED_VALUE).description("Observed value of the original series for the period where new data was added").outputClass(double[].class).build())
        .item(AtomicDictionary.Item.builder().name(SERIES_NEWS_T).description("News of the standardized series for the period where new data was added").outputClass(double[].class).build())
        .item(AtomicDictionary.Item.builder().name(SERIES_NEWS).description("News of the original series for the period where new data was added").outputClass(double[].class).build())
        .item(AtomicDictionary.Item.builder().name(SERIES_WEIGHTS_T).description("Weights of the news of each standardized series on forecasts of a target series").outputClass(Matrix.class).build())
        .item(AtomicDictionary.Item.builder().name(SERIES_WEIGHTS).description("Weights of the news of each original series on forecasts of a target series").outputClass(Matrix.class).build())
        .item(AtomicDictionary.Item.builder().name(SERIES_IMPACTS_T).description("Total impacts of the news of each standardized series on forecasts of a target series").outputClass(Matrix.class).build())
        .item(AtomicDictionary.Item.builder().name(SERIES_IMPACTS).description("Total impacts of the news of each original series on forecasts of a target series").outputClass(Matrix.class).build())
        .item(AtomicDictionary.Item.builder().name(FORECASTS_PERIODS).description("Forecasts periods").outputClass(String[].class).build())
        .item(AtomicDictionary.Item.builder().name(OLD_FORECASTS_T).description("Previous forecasts of the standardized target series").outputClass(double[].class).build())
        .item(AtomicDictionary.Item.builder().name(OLD_FORECASTS).description("Previous forecasts of the original target series").outputClass(double[].class).build())
        .item(AtomicDictionary.Item.builder().name(REVISED_FORECASTS_T).description("Forecasts of the standardized target series based on revised data").outputClass(double[].class).build())
        .item(AtomicDictionary.Item.builder().name(REVISED_FORECASTS).description("Forecasts of the original target series based on revised data").outputClass(double[].class).build())
        .item(AtomicDictionary.Item.builder().name(NEW_FORECASTS_T).description("New forecasts of the standardized target series").outputClass(double[].class).build())
        .item(AtomicDictionary.Item.builder().name(NEW_FORECASTS).description("New forecasts of the original target series").outputClass(double[].class).build())
        .build();      
}
