# Copyright (c) 2017. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import warnings

class DataFrameLoader(object):
    """
    Wraps a `DataFrame` with some information on how to join it.

    Parameters
    __________
    name : str
        The name of the dataframe, to easily reference it.
    load_dataframe : function
        A function that returns the `DataFrame` object.
    join_on_right : str
        The column of the `DataFrame` to join on (i.e. the patient
        ID column name).
    join_on_left: str
        The corresponding column in the `cohorts.Patient.additional_data`
	to join on, if not the `id` (which is the default).
    """
    def __init__(self,
                 name,
                 load_dataframe,
                 join_on=None,
                 join_on_right="patient_id",
                 join_on_left="patient_id"):
        self.name = name
        self.load_dataframe = load_dataframe
        if join_on is not None:
            warnings.warn("`join_on` parameter is deprecated. Please use `join_on_right` instead.", DeprecationWarning)
            self.join_on_right = join_on
        else:
            self.join_on_right = join_on_right
        self.join_on_left = join_on_left
