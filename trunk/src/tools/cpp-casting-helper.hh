/*
    Copyright(C) 2006 Frans Englich <frans.englich@telia.com>

    This file is part of the KDE project

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Library General Public
    License as published by the Free Software Foundation; either
    version 2 of the License, or(at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Library General Public License for more details.

    You should have received a copy of the GNU Library General Public License
    along with this library; see the file COPYING.LIB.  If not, write to
    the Free Software Foundation, Inc., 51 Franklin Street - Fifth Floor,
    Boston, MA 02110-1301, USA.
*/

#ifndef Patternist_CppCastingHelper_H
#define Patternist_CppCastingHelper_H

namespace Patternist
{
    /**
     * @short Provides convenience methods for performing static casts between C++ classes.
     *
     * In Patternist, it is very common to do up-casts from Expression or Item, which typically
     * involves writing messy code. Such an old-way cast looks like this:
     *
     * @code
     * static_cast<const MyClass *>(myInstance.get())->myClassMember()
     * @endcode
     * 
     * CppCastingHelper provides the convenience method as() for this, which is functionally
     * equivalent to the above code, but simpler:
     * 
     * @code
     * myInstance->as<MyClass>()->myClassMember()
     * @endcode
     *
     * The as() function performs a static cast.
     * 
     * By using CppCastingHelper, this is achieved:
     *
     * - Const correctness is automatically taken care of
     * - Less code to write
     * - When compiling in debug mode, the as() functions uses a @c dynamic_cast to verify that the
     *   static casts are properly done, such that sensible error messages are given when the casts
     *   are invalid. It also traps invalid casts which nevertheless happen to work on a particular
     *   platform/compiler/hardware architecture.
     *
     * CppCastingHelper is a template class where the TSubClass parameter must be the class
     * inheriting CppCastingHelper. See Item or Expression for demonstration.
     *
     * @author Frans Englich <frans.englich@telia.com>
     */
    template<typename TSubClass>
    class CppCastingHelper
    {
    public:

        /**
         * Casts this instance to:
         *
         * @code
         * const TCastTarget *
         * @endcode
         *
         * and returns the result.
         *
         * When compiled in debug mode, this function perform a @c dynamic_cast, in order to
         * produce an informative message.
         */
        template<typename TCastTarget>
        inline const TCastTarget *as() const
        {
            Q_ASSERT_X(dynamic_cast<const TCastTarget *>(static_cast<const TSubClass *>(this)),
                       "CppCastingHelper::as() const",
                       "The cast is invalid. This class does not inherit the cast target.");
            return static_cast<const TCastTarget *>(static_cast<const TSubClass *>(this));
        }

        /**
         * Casts this instance to:
         *
         * @code
         * TCastTarget *
         * @endcode
         *
         * and returns the result.
         *
         * When compiled in debug mode, this function perform a @c dynamic_cast, in order to
         * produce an informative message.
         */
        template<typename TCastTarget>
        inline TCastTarget *as()
        {
            Q_ASSERT_X(dynamic_cast<TCastTarget *>(static_cast<TSubClass *>(this)),
                       "CppCastingHelper::as()",
                       "The cast is invalid. This class does not inherit the cast target.");
            return static_cast<TCastTarget *>(static_cast<TSubClass *>(this));
        }

    protected:
        /**
         * This constructor is protected because this class must be sub-classed.
         */
        inline CppCastingHelper() {}
    };
}

#endif
// vim: et:ts=4:sw=4:sts=4
