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
import jdplus.toolkit.base.api.information.InformationExtractor;

module jdplus.dfm.base.core {
    
    requires static lombok;
    requires static nbbrd.design;
    requires static nbbrd.service;
    requires static org.checkerframework.checker.qual;

    requires transitive jdplus.toolkit.base.core;
    requires transitive jdplus.dfm.base.api;

    exports jdplus.dfm.base.core;
    exports jdplus.dfm.base.core.var;
    exports jdplus.dfm.base.core.varma;
    
    provides InformationExtractor with
        jdplus.dfm.base.core.extractors.DynamicFactorModelExtractor,
        jdplus.dfm.base.core.extractors.DfmEstimatesExtractor,
        jdplus.dfm.base.core.extractors.DfmResultsExtractor,
        jdplus.dfm.base.core.extractors.DfmResultsNewsExtractor;
}
