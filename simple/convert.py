import sys
from dump import dump
from ensight import ensight
d = dump("dump_t1");
d.map(1,"id",2,"type",3,"x",4,"y",5,"z");
e = ensight(d);
e.one("PMBt1")
