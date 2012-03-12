/**
\page igeslib GoTools Igeslib

The \beginlink \link Go::IGESconverter IGES converter \endlink
read an IGES file and represents its entites in the 
internal data structure of GoTools. It can also write a model represented in
GoTools to an IGES file or convert between an IGES file and the 
\beginlink \link streamable_doc internal file format \endlink of GoTools.

GoTools represent only geometric entities. Thus, IGES entities like annotation,
structure, property, associativity, view, drawing 
and figure will be neglected. Neither are constructive solid geometry or
finite element modelling entites handled. If such entities exist in a file
read by the IGES converter, warning messages will be issued.

The topological entities specified in IGES 5.3 is not handled by the current
version of the IGES converter. Thus, the entities vertex, edge, edge list,
loop, face and shell is not handled. However, the geometric entities 
corresponding to these topological entities will be read. Colour information
is read.

The content of an IGES file is transferred to the application as a vector
of \beginlink \link Go::GeomObject GeomObjects\endlink.
By checking the type of each 
object and acting thereafter, the model can be stored and handled in 
the GoTools environment.

To write an IGES file, the file entities are added one by one to the IGES
convertor using the function addGeom which takes a GeomObject as parameter.
The actual file is written by the command writeIGES.

*/
