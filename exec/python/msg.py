# @file:    msg.py
# @purpose: message logging
# @author:  Doug Schouten <dschoute@sfu.ca>

import sys

#_______________________________________________________________________________
class msglog:
    levels = { 'debug'   : 0,
               'info'    : 1,
               'warning' : 2,
               'error'   : 3,
               'fatal'   : 4 }
    def __init__( self, name, level = 'debug', useColor = True ):
        self._level = self.levels[level.lower()]
        self._levelname = level
        self._name = name
        self._printonly = True
        self._usecolor = useColor
        
    def _print( self, msg, level ):
        print '%s from %s %s'%( level.upper(), self._name, msg )
                              
    def setLevel( self, level ):
        self._level = self.levels[ level.lower() ]
        self._levelname = level.lower()

    def setPrintOnly( self, v = True ):
        self._printonly = v

    def toString( self, args ):
        return ' '.join( map( str, args ) )

    def debug( self, *msg ):
        if self._level <= self.levels['debug']:
            if self._usecolor:
                self._print( '\033[1;30m'+self.toString(msg)+'\033[0m', 'debug' )
            else:
                self._print( self.toString(msg), 'debug' )
        pass
    
    def info( self, *msg ):
        if self._level <= self.levels['info']:
            if self._usecolor:
                self._print( '\033[1;32m'+self.toString(msg)+'\033[0m', 'info' )
            else:
                self._print( self.toString(msg), 'info' )
        pass
    
    def warning( self, *msg ):
        if self._level <= self.levels['warning']:
            if self._usecolor:
                self._print( '\033[1;33m'+self.toString(msg)+'\033[0m', 'warning' )
            else:
                self._print( self.toString(msg), 'warning' )
        pass

    def error( self, *msg ):
        if self._level <= self.levels['error']:
            if self._usecolor:
                self._print( '\033[1;31m'+self.toString(msg)+'\033[0m', 'error' )
            else:
                self._print( self.toString(msg), 'error' )
        if not self._printonly:
            raise RuntimeError
        pass

    def fatal( self, *msg ):
        if self._level <= self.levels['fatal']:
            if self._usecolor:
                self._print( '\033[1;48m'+self.toString(msg)+'\033[1;48m', 'fatal' )
            else:
                self._print( self.toString(msg), 'fatal' )
        if not self._printonly:
            sys.exit( 1 )
        pass
    
