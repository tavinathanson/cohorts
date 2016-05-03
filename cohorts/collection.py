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

class Collection(object):
    def __init__(self, elements):
        self.elements = elements

    def short_string(self):
        """
        Compact string representation which doesn't print any of the
        collection elements.
        """
        file_str = ""
        return "<%s with %d elements>" % (
            self.__class__.__name__,
            len(self))

    def to_string(self, limit=None):
        """
        Create a string representation of this collection, showing up to
        `limit` items.
        """
        header = self.short_string()
        if len(self) == 0:
            return header
        contents = ""
        element_lines = [
            "  -- %s" % (element,)
            for element in self.elements[:limit]
        ]
        contents = "\n".join(element_lines)

        if limit is not None and len(self.elements) > limit:
            contents += "\n  ... and %d more" % (len(self) - limit)
        return "%s\n%s" % (header, contents)

    def __str__(self):
        return self.to_string(limit=50)

    def __repr__(self):
        return str(self)

    def __len__(self):
        return len(self.elements)

    def __iter__(self):
        return iter(self.elements)

    def __getitem__(self, idx):
        return self.elements[idx]
