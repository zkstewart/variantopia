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
from contextlib import contextmanager

@contextmanager
def annotarium_import_context(annotariumDir):
    '''
    Context manager for importing from annotarium. annotarium's modules import
    each other by bare name (e.g., 'from validation import ...'), and some of
    those names are shared with variantopia's modules; while importing, those
    variantopia modules are set aside and annotarium's modules directory is
    searched first, so that annotarium gets its own modules.
    '''
    modulesDir = os.path.join(annotariumDir, "modules")
    clashing = {
        name: module for name, module in sys.modules.items()
        if "." not in name and os.path.isfile(os.path.join(modulesDir, name + ".py"))
    }
    for name in clashing:
        del sys.modules[name]
    sys.path.insert(0, modulesDir)
    sys.path.append(os.path.dirname(annotariumDir))
    try:
        yield
    finally:
        sys.path.remove(modulesDir)
        sys.modules.update(clashing)

def import_annotarium_gff3(annotariumDir):
    try:
        with annotarium_import_context(annotariumDir):
            from annotarium import GFF3Feature, GFF3Tarium
        return GFF3Feature, GFF3Tarium
    except (ModuleNotFoundError, ImportError):
        raise ModuleNotFoundError(f"Could not import GFF3 classes from '{annotariumDir}'")

def import_annotarium_domains(annotariumDir):
    try:
        with annotarium_import_context(annotariumDir):
            from annotarium import Domains, OverlapResolver
        return Domains, OverlapResolver
    except (ModuleNotFoundError, ImportError):
        raise ModuleNotFoundError(f"Could not import Domains classes from '{annotariumDir}'")

def import_annotarium_annotable(annotariumDir):
    try:
        with annotarium_import_context(annotariumDir):
            from annotarium import Annotable
        return Annotable
    except (ModuleNotFoundError, ImportError):
        raise ModuleNotFoundError(f"Could not import the Annotable class from '{annotariumDir}'")
