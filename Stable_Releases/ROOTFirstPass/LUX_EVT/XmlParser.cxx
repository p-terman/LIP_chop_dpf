//
//  Taken and modified from $ROOTSYS/tutorials/xml/xmlreadfile.C
//

#include <stdlib.h>
#include <iostream>
#include <string>

#include "TXMLEngine.h"
#include "XmlParser.h"

XmlParser::XmlParser(const char* buffer)
{
    // First create engine
    TXMLEngine* xml = new TXMLEngine;
   
    // Now try to parse xml file
    // Only file with restricted xml syntax are supported
    XMLDocPointer_t xmldoc = xml->ParseString( buffer );
    if (xmldoc==0) {
        delete xml;
        std::cout << "Error in XML Header: cannot Parse! " << std::endl;
        return;  
    }

    // take access to main node   
    XMLNodePointer_t mainnode = xml->DocGetRootElement(xmldoc);
   
    // display recursively all nodes and subnodes
    DisplayNode(xml, mainnode);

    // added to read all <settings> tags, not just the first and its children
    XMLNodePointer_t nextNode = xml->GetNext(mainnode);
    while (nextNode!=0) {
        DisplayNode(xml, nextNode);
        nextNode = xml->GetNext(nextNode);
    }    
    // Release memory before exit
    xml->FreeDoc(xmldoc);
    delete xml;
}


void XmlParser::DisplayNode( TXMLEngine* xml, XMLNodePointer_t node ) 
{
    std::string nodeName;
    std::string parentName;
    nodeName = xml->GetNodeName(node);

    if( nodeName.find("daq_version") < std::string::npos ) {
        globalDaq = atof( xml->GetNodeContent(node) ) ;
    }

    if( nodeName.find("event_builder_version") < std::string::npos ) { 
        parentName = xml->GetNodeName(xml->GetParent(node));
        if( parentName.find("evt_settings") < std::string::npos )
            evtSettingsEventBuilder = atof( xml->GetNodeContent(node));
        else
            globalEventBuilder = atof( xml->GetNodeContent(node) ) ;
   }
     
    // display namespace  for <node:namespace>content</node:namespace>
    // xml->GetNS(node);
    // if (ns!=0)
   
    // display attributes  for <node>attr1="val1" attr2="val2"</node>
    // XMLAttrPointer_t attr = xml->GetFirstAttr(node);
    // while (attr!=0) {
    //   attr = xml->GetNextAttr(attr);  
    // }
   
    // display all child nodes   
    XMLNodePointer_t child = xml->GetChild(node);
    while (child!=0) {
        DisplayNode(xml, child); 
        child = xml->GetNext(child);
    }
}
