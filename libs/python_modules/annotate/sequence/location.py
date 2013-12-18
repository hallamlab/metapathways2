"""This module defines classes for working with locations and coordinates."""

import re
class Location(object):
    """A location containing at least one pair of :class:`Coordinates`
    *coordinates*."""

    def __init__(self, coordinates):
        self.coordinates = coordinates

    # XXX FIXME This is a bug: join() is connecting complement() sections together: 
    def __str__(self):
        if len(self.coordinates) > 1:
            return 'join(%s)' % ','.join([str(c) for c in self.coordinates])
        return str(self.coordinates[0])

    def __get_begin(self):
        return min(int(c.begin) for c in self.coordinates)
    #: The minimum start coordinate of this location feature
    begin = property(__get_begin)
    start = property(__get_begin)

    def __get_end(self):
        return max(int(c.end) for c in self.coordinates)
    #: The maximum end coordinate of this location feature
    end = property(__get_end)

    def is_complement_strand(self):
        return self.begin > self.end
    #: Whether the location is a reference to the complement strand
    complement = property(is_complement_strand)



class Coordinates(object):
    """A pair of coordinates starting at the :class:`Coordinate` *begin* and
    ending at the :class:`Coordinate` *end*."""

    def __init__(self, begin, end):
        self.begin = begin
        self.end = end
    
    def __str__(self):
        if self.begin > self.end:
            return 'complement(%s..%s)' % (self.end, self.begin)
        return '%s..%s' % (self.begin, self.end)

    def __repr__(self):
        return '<Coordinates %s..%s>' % (repr(self.begin), repr(self.end))



numberPATTERN = re.compile(r"(\d+)")
class Coordinate(object):
    """A coordinate at *coordinate*."""
    def __init__(self, coordinate):

        m = numberPATTERN.search(coordinate)
        if m:
           extracted_coordinate = m.group(1)
        else:
           extracted_coordinate = 0

        self.coordinate = int( extracted_coordinate )

    def __int__(self):
        return self.coordinate

    def __cmp__(self, other):
        return self.coordinate - other.coordinate

    def __str__(self):
        return '%d' % self.coordinate

    def __repr__(self):
        return '<Coordinate %s>' % self.coordinate



class FuzzyCoordinate(Coordinate):
    """An imprecise coordinate centred around *coordinate*, with a range
    extended by *bound*."""
    def __init__(self, coordinate, bound='<'):
        m = numberPATTERN.search(coordinate)
        if m:
           extracted_coordinate = m.group(1)
        else:
           extracted_coordinate = 0

        self.coordinate = int( extracted_coordinate )
        self.bound = bound
 

    def __str__(self):
        return '%s%d' % (self.bound, self.coordinate)

    def __repr__(self):
        return '<FuzzyCoordinate %s%s>' % (self.bound, self.coordinate)
