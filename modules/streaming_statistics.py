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

from statistics import median, StatisticsError

class Remedian:
    '''
    This function is adapted from code originally found in https://ttv1.github.io/remedian.html,
    with the _medianPrim function coming from https://github.com/emilianavt/OpenSeeFace/blob/master/remedian.py
    as it was not originally described in the ttv1 implementation (despite the placeholder being there).
    
    As such, this function should come with the licence originally associated with the ttv1 implementation,
    and it is embedded here as:
    
    Copyright © 2016,2017 Tim Menzies tim@menzies.us, MIT license v2.
    
    Permission is hereby granted, free of charge, to any person obtaining a copy of this software
    and associated documentation files (the "Software"), to deal in the Software without restriction,
    including without limitation the rights to use, copy, modify, merge, publish, distribute,
    sublicense, and/or sell copies of the Software, and to permit persons to whom the Software
    is furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in all copies or substantial
    portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT
    NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
    IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
    WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
    SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
    '''
    def __init__(self, inits=[], k=64):
        self.all = []
        self.k = k
        self.more = None
        self._median = None
        self._count = 0
        [self + x for x in inits]
    
    def __add__(self, x):
        '''
        This method allows an object instance to be used as 'remedianObj + 5' to add a
        new number to it. Do not use 'remedianObj += 5' as that will return None.
        '''
        self._count += 1
        self._median = None
        self.all.append(x)
        if len(self.all) == self.k:
            self.more = self.more or Remedian(k=self.k)
            self.more + self._medianPrim(self.all)
            self.all = []  # reset
    
    @property
    def median(self):
        if self._count == 0:
            return None
        
        return self.more.median if self.more else self._medianPrim(self.all)
    
    def _medianPrim(self, all):
        if self._median == None:
            self._median = median(all)
        return self._median
    
    def __repr__(self):
        return "<Remedian object;k={0};median={1}>".format(
            self.k,
            self.median
        )

class OnlineMean:
    def __init__(self, inits=[]):
        self.mean = None
        self._count = 0
        [self + x for x in inits]
    
    def __add__(self, x):
        '''
        This method allows an object instance to be used as 'meanObj + 5' to add a
        new number to it. Do not use 'meanObj += 5' as that will return None.
        '''
        self._count += 1
        if self.mean == None:
            self.mean = 0
        self.mean += (x - self.mean) / self._count
    
    def __repr__(self):
        return "<OnlineMean object;mean={0}>".format(
            self.mean
        )

def merge_online_means(inputs):
    totalCount = 0
    summedMeans = 0
    for onlineMeanObj in inputs:
        if onlineMeanObj.mean != None:
            summedMeans += onlineMeanObj.mean * onlineMeanObj._count
            totalCount += onlineMeanObj._count
    
    if totalCount == 0:
        return None
    else:
        return summedMeans / totalCount

def merge_remedians(inputs):
    medians = []
    for remedianObj in inputs:
        if remedianObj.median != None:
            medians.append(remedianObj.median)
    
    if len(medians) == 0:
        return None
    else:
        return median(medians)
