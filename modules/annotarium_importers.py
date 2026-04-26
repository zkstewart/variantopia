# Copyright (C) 2026 Zachary Kenneth Stewart

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

import os, sys

def import_annotarium_gff3(annotariumDir):
    try:
        sys.path.append(os.path.dirname(annotariumDir))
        from annotarium import GFF3Feature, GFF3Tarium
        return GFF3Feature, GFF3Tarium
    except ModuleNotFoundError:
        raise ModuleNotFoundError(f"Could not import GFF3 classes from '{args.annotariumDir}'")

def import_annotarium_domains(annotariumDir):
    try:
        sys.path.append(os.path.dirname(annotariumDir))
        from annotarium import Domains, OverlapResolver
        return Domains, OverlapResolver
    except ModuleNotFoundError:
        raise ModuleNotFoundError(f"Could not import Domains classes from '{args.annotariumDir}'")
