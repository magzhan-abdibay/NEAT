cmake_minimum_required(VERSION 2.4)

PROJECT( rtneat )

SET(rtneat_sources 
  experiments.cpp
        Gene.cpp
        Genome.cpp
        Innovation.cpp
        Link.cpp
        Neat.cpp
  neatmain.cpp
        Network.cpp
        NNode.cpp
        Organism.cpp
        Population.cpp
        Species.cpp
        Trait.cpp
  CartPole.cpp)
#FILE(GLOB rtneat_sources *.cpp)

ADD_EXECUTABLE( rtneat ${rtneat_sources} )

