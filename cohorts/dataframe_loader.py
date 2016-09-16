# Copyright (c) 2016. Mount Sinai School of Medicine
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

from __future__ import print_function

class DataFrameLoader(object):
    """
    Wraps a `DataFrame` with some information on how to join it.

    Parameters
    __________
    name : str
        The name of the dataframe, to easily reference it.
    load_dataframe : function
        A function that returns the `DataFrame` object.
    join_on : str
        The column of the `DataFrame` to join on (i.e. the patient
        ID column name).
    """
    def __init__(self,
                 name,
                 load_dataframe,
                 join_on="patient_id"):
        self.name = name
        self.load_dataframe = load_dataframe
        self.join_on = join_on
