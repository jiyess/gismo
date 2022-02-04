/** @file gsXmlCollection.h

    @brief Provides declaration of input/output XML utilities struct.

    This file is part of the G+Smo library. 

    This Source Code Form is subject to the terms of the Mozilla Public
    License, v. 2.0. If a copy of the MPL was not distributed with this
    file, You can obtain one at http://mozilla.org/MPL/2.0/.
    
    Author(s): A. Mantzaflaris
*/

#pragma once

#include <gsCore/gsLinearAlgebra.h>

#include <iostream>
#include <string>

#include <gsIO/gsXml.h>
#include <gsIO/gsFileData.h>

namespace gismo {

class gsXmlCollection
{
public:
    typedef std::string String;
public:

    /// Constructor using a filename.
    gsXmlCollection(String const& fn)
    : m_fn(fn)
    {
        // mfile <<"<?xml version=\"1.0\"?>\n";
        // mfile <<"<VTKFile type=\"Collection\" version=\"0.1\">";
        // mfile <<"<Collection>\n";

    }

    /// Adds a part in the collection, with complete filename (including extension) \a fn
    void addFile(String const & fn)
    {
        // GISMO_ASSERT(fn.find_last_of(".") != String::npos, "File without extension");
        // GISMO_ASSERT(counter!=-1, "Error: collection has been already saved." );
        // mfile << "<DataSet part=\""<<counter++<<"\" file=\""<<fn<<"\"/>\n";
        m_fd.addString(fn);
    }

    /// Adds a part in the collection, with complete filename (including extension) \a fn
    void addFile(String const & label, String const & fn)
    {
        // GISMO_ASSERT(fn.find_last_of(".") != String::npos, "File without extension");
        // GISMO_ASSERT(counter!=-1, "Error: collection has been already saved." );
        // mfile << "<DataSet part=\""<<counter++<<"\" file=\""<<fn<<"\"/>\n";
        m_fd.addString(fn, label);
    }

    /// Finalizes the collection by closing the XML tags, always call
    /// this function (once) when you finish adding files
    void save()
    {
        // GISMO_ASSERT(counter!=-1, "Error: gsParaviewCollection::save() already called." );
        // mfile <<"</Collection>\n";
        // mfile <<"</VTKFile>\n";

        // mfn.append(".pvd");
        // std::ofstream f( mfn.c_str() );
        // GISMO_ASSERT(f.is_open(), "Error creating "<< mfn );
        // f << mfile.rdbuf();
        // f.close();
        // mfile.str("");
        // counter = -1;
        m_fd.save(m_fn.str());
    }

    /// Finalizes the collection by closing the XML tags, always call
    /// this function (once) when you finish adding files
    void load(String const & fn)
    {
        m_fd.clear();
        m_fd.read(fn);
    }

    /// Get Functions
public:

    template<class Object>
    inline void getId( const int & id, Object& result)  const
    {
        gsFileData<> file(m_fd.getString(id));
        file.getFirst(result);
    }

    template<class Object>
    inline void getLabel( const std::string & label, Object& result)  const
    {
        //gsFileData<> file(m_fd.getString(id));
        //gsInfo << m_fd.getString(id) << "\n";
        //file.getFirst(result);
        gsFileData<> file(m_fd.getString(label));
        file.getFirst(result);
    }

    gsMatrix<real_t> getMatrix( const std::string & label)  const
    {
        gsMatrix<real_t> result;
        gsFileData<real_t> file(m_fd.getString(label));
        file.getFirst(result);
        return result;
    }

    gsMatrix<real_t> getMatrix(  const int & id )  const
    {
        gsMatrix<real_t> result;
        gsFileData<real_t> file(m_fd.getString(id));
        file.getFirst(result);
        return result;
    }

    std::string getString( index_t  id)  const
    {
        return m_fd.getString(id);
    }

    /// \brief Reports the number of objects which are held in the file data
    int numData() const { return m_fd.numData();}


private:
    /// Pointer to char stream
    std::stringstream m_fn;
    gsFileData<> m_fd;

};

#ifdef GISMO_BUILD_PYBIND11

  /**
   * @brief Initializes the Python wrapper for the class: gsXmlCollection
   */
  void pybind11_init_gsXmlCollection(pybind11::module &m);

#endif // GISMO_BUILD_PYBIND11

}// end namespace gismo

#ifndef GISMO_BUILD_LIB
#include GISMO_HPP_HEADER(gsXmlCollection.hpp)
#endif