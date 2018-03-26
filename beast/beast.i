%module beast
%{
#include "config.h"
#include "stars.h"
#include "constellations.h"
#include "beast.h"
%}
//newobject gives python control of these objects
%newobject star_db::copy;
%newobject star_db::copy_n_brightest;
%newobject star_db::search;
%newobject star_db::operator-;
%newobject star_db::operator&;
%newobject star_query::from_kdmask;
%newobject star_query::from_kdresults;
%newobject match_result::from_match;
%include "config.h"
%include "stars.h"
%include "constellations.h"
%include "beast.h"
