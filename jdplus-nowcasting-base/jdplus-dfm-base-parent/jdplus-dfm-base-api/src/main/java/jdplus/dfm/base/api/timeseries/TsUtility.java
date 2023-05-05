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
package jdplus.dfm.base.api.timeseries;

import java.time.LocalDate;
import jdplus.toolkit.base.api.timeseries.TsData;
import jdplus.toolkit.base.api.timeseries.TsPeriod;
import jdplus.toolkit.base.api.timeseries.TsUnit;

/**
 *
 * @author palatej
 */
@lombok.experimental.UtilityClass
public class TsUtility {

    public TsPeriod lastPeriod(TsPeriod p, int annualFrequency) {
        int pfreq = p.annualFrequency();
        if (pfreq == annualFrequency) {
            return p;
        }
        TsUnit unit = TsUnit.ofAnnualFrequency(annualFrequency);
        LocalDate last = p.end().toLocalDate().minusDays(1);
        return TsPeriod.of(unit, last);
    }

    public TsPeriod firstPeriod(TsPeriod p, int annualFrequency) {
        int pfreq = p.annualFrequency();
        if (pfreq == annualFrequency) {
            return p;
        }
        TsUnit unit = TsUnit.ofAnnualFrequency(annualFrequency);
        return TsPeriod.of(unit, p.start());
    }

    public boolean isMissing(TsData s, int idx) {
        return Double.isFinite(s.getValue(idx));
    }

    public LocalDate lastDay(TsPeriod p) {
        return p.end().toLocalDate().minusDays(1);
    }

    public TsData extend(final TsData data, final LocalDate lastday) {
        TsPeriod s = TsPeriod.of(data.getTsUnit(), lastday);
        if (!lastday.equals(lastDay(s))) {
            s = s.previous();
        }
        int n = data.getStart().until(s) + 1;
        return data.extend(0, n - data.length());
    }

}
