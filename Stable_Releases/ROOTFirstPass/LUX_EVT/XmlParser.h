#ifndef XML_H
#define XML_H 1

#include "TXMLEngine.h"

class XmlParser {
    public:
        XmlParser() {};
        XmlParser(const char* );
        ~XmlParser() {};
        
    public:
        void DisplayNode(TXMLEngine* , XMLNodePointer_t );
        double GetEvtSettingsEventBuilderVersion() { return evtSettingsEventBuilder; };
        double GetGlobalEventBuilderVersion() { return globalEventBuilder; };
        double GetGlobalDAQVersion() { return globalDaq; };

    private:
        double evtSettingsEventBuilder ;
        double globalEventBuilder ;
        double globalDaq;

};

#endif
