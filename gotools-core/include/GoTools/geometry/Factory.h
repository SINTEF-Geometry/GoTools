//===========================================================================
//                                                                           
// File: Factory.h                                                           
//                                                                           
// Created: Mon Nov  6 10:00:43 2000                                         
//                                                                           
// Author: Atgeirr F Rasmussen <atgeirr@sintef.no>
//                                                                           
// Revision: $Id: Factory.h,v 1.14 2009-05-13 07:30:48 vsk Exp $
//                                                                           
// Description:
//                                                                           
//===========================================================================

#ifndef _FACTORY_H
#define _FACTORY_H

#include <map>
#include "GoTools/geometry/GeomObject.h"
#include <memory>
#include "GoTools/utils/config.h"

namespace Go
{
    /// Abstract base class for Creators.  These are policy classes for generating
    /// objects derived from GeomObject.  This class is only intended for internal use 
    /// by the \ref Factory class.  The user should not have to worry about it.
    class Creator
    {
    public:
	virtual GeomObject* create() = 0;
	virtual ~Creator() {}
    };

    /// This is the concrete Creator class for generating GeomObject-derived classes.
    /// The class is only intended for internal use by the \ref Factory class.  The user
    /// should not have to worry about it.
    template <class T>
    class ConcreteCreator : public Creator
    {
    public:
	ConcreteCreator() {}
	GeomObject* create()
	{
	    return new T;
	}
    };
    
    /** This is the Factory for creating instances of GeomObjects of type
     *  requested by the user.  There will be one, global Factory object,
     *  called Go::global_factory_.  However, the user will not have to deal
     *  with this directly.  The user need only interact with the Factory 
     *  using the static member function 'createObject(ClassType)' and the
     *  template Register() function (or Registrator class, if the compiler
     *  misbehaves. In practice, this means calling 'GoTools::init()' early
     *  in the code.
     */
    class Factory
    {
    public:
	~Factory()
	{
	    std::map<ClassType, Creator*>::iterator it;
	    for (it = themap_.begin(); it != themap_.end(); ++it)
		delete it->second;
	}

	/// Makes a new GeomObject instance of the specified ClassType and returns
	/// a pointer to it.  The user assumes ownership over the created object.
	/// \param class_type the class type of the object that the user wants to have
	///                   constructed.
	static GeomObject* createObject(ClassType class_type)
	{
	    return globalFactory()->doCreateObject(class_type);
	}

	/// Register a ClassType with the Factory.  This amounts to provide the Factory
	/// with the Creator object that is used to generate a new GeomObject of type
	/// ClassType.  Usually, the user would not want to call this function directly,
	/// but through the (template) function \ref Register() or by using a
	/// \ref Registrator
	/// \param class_type the ClassType for which we want the Factory to associate
	///                   a particular Creator.
	/// \param c pointer to the Creator object to be used by the Factory when creating
	///          objects of type ClassType.
	void registerClass(ClassType class_type, Creator* c)
	{
	    themap_[class_type] = c;
	}

	/// This function returns a pointer to the unique, global Factory.
	/// \return pointer to global Factory.
	static Factory* globalFactory()
	{
	    static Factory factory_singleton;
	    return &factory_singleton;
// 	    if (global_factory_.get() == 0)
// 		global_factory_ = std::shared_ptr<Factory>(new Factory);
// 	    return global_factory_.get();
	}
    private:
	// private constructor.  User should not need to explicitly construct 
	// Factory.
	Factory()
	{
	}

	GeomObject* doCreateObject(ClassType class_type)
	{
	    std::map<ClassType, Creator*>::iterator it;
	    it = themap_.find(class_type);
	    if (it==themap_.end()) {
		THROW("Class type number " << class_type
		      << " is not registered.");
	    }
	    return it->second->create();
	}
	//static GO_API std::shared_ptr<Factory> global_factory_;
	std::map<ClassType, Creator*> themap_;
    };

    /// This function is used to register a class derived from GeomObject
    /// with the global Factory.  By using this function rather than
    /// \ref Factory::registerClass(), the user does not have to worry about the
    /// details of the Creator class.  
    /// To register a class \code DerivedClass \endcode, it should be 
    /// sufficient to run: \code Register<DerivedClass>() \endcode .
    template <class T>
    void Register()
    {
	Factory* f = Factory::globalFactory();
	ConcreteCreator<T>* c = new ConcreteCreator<T>;
	f->registerClass(T::classType(), c);
	delete c;
    }

    /** Work around for compilation problems.
     *
     */

    // @afr: I have no idea why, but VS6 refuses to do 

    /// On some compilators (ie., VS6), the \ref Register() function cannot
    /// be used directly.  To work around this, we wrap it in the constructor
    /// of a dummy class, which we call registrator.  With other words, in 
    /// order to register the class 'DerivedClass' (supposedly inheriting from
    /// GeomObject), we type: \code Registrator<DerivedClass> r \endcode.
    template <class T>
    class Registrator
    {
    public:
	Registrator()
	{
	    Factory* f = Factory::globalFactory();
	    ConcreteCreator<T>* c = new ConcreteCreator<T>;
	    f->registerClass(T::classType(), c);
	}
    };


} // namespace Go

#endif // _FACTORY_H

