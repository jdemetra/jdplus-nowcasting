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

import jdplus.dfm.base.core.IDfmMeasurement;
import jdplus.toolkit.base.core.data.DataBlock;


/**
 * Z = 1 [0 ... 0]
 *
 * @author Jean Palate
 */
public class LevelMeasurement implements IDfmMeasurement {

    public static final LevelMeasurement ML = new LevelMeasurement();

    private LevelMeasurement() {
    }

    @Override
    public int getLength() {
        return 1;
    }

    @Override
    public void fill(DataBlock z) {
        z.set(0, 1);
    }

    @Override
    public double dot(DataBlock x) {
        return x.get(0);
    }

}
