/* Software License Agreement (BSD License)
 *
 * Copyright (c) 2011, Willow Garage, Inc.
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *
 *  * Redistributions of source code must retain the above copyright
 *    notice, this list of conditions and the following disclaimer.
 *  * Redistributions in binary form must reproduce the above
 *    copyright notice, this list of conditions and the following
 *    disclaimer in the documentation and/or other materials provided
 *    with the distribution.
 *  * Neither the name of Willow Garage, Inc. nor the names of its
 *    contributors may be used to endorse or promote products derived
 *    from this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
 * BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES;
 * LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
 * CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
 * ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * Auto-generated by genmsg_cpp from file /home/justin/charm/charm/trunk/ROS/catkin_ws/src/wam_msgs/msg/RTOrtnVel.msg
 *
 */


#ifndef WAM_MSGS_MESSAGE_RTORTNVEL_H
#define WAM_MSGS_MESSAGE_RTORTNVEL_H


#include <string>
#include <vector>
#include <map>

#include <ros/types.h>
#include <ros/serialization.h>
#include <ros/builtin_message_traits.h>
#include <ros/message_operations.h>


namespace wam_msgs
{
template <class ContainerAllocator>
struct RTOrtnVel_
{
  typedef RTOrtnVel_<ContainerAllocator> Type;

  RTOrtnVel_()
    : angular()
    , magnitude(0.0)  {
      angular.assign(0.0);
  }
  RTOrtnVel_(const ContainerAllocator& _alloc)
    : angular()
    , magnitude(0.0)  {
      angular.assign(0.0);
  }



   typedef boost::array<float, 3>  _angular_type;
  _angular_type angular;

   typedef float _magnitude_type;
  _magnitude_type magnitude;




  typedef boost::shared_ptr< ::wam_msgs::RTOrtnVel_<ContainerAllocator> > Ptr;
  typedef boost::shared_ptr< ::wam_msgs::RTOrtnVel_<ContainerAllocator> const> ConstPtr;
  boost::shared_ptr<std::map<std::string, std::string> > __connection_header;

}; // struct RTOrtnVel_

typedef ::wam_msgs::RTOrtnVel_<std::allocator<void> > RTOrtnVel;

typedef boost::shared_ptr< ::wam_msgs::RTOrtnVel > RTOrtnVelPtr;
typedef boost::shared_ptr< ::wam_msgs::RTOrtnVel const> RTOrtnVelConstPtr;

// constants requiring out of line definition



template<typename ContainerAllocator>
std::ostream& operator<<(std::ostream& s, const ::wam_msgs::RTOrtnVel_<ContainerAllocator> & v)
{
ros::message_operations::Printer< ::wam_msgs::RTOrtnVel_<ContainerAllocator> >::stream(s, "", v);
return s;
}

} // namespace wam_msgs

namespace ros
{
namespace message_traits
{



// BOOLTRAITS {'IsFixedSize': True, 'IsMessage': True, 'HasHeader': False}
// {'wam_msgs': ['/home/justin/charm/charm/trunk/ROS/catkin_ws/src/wam_msgs/msg']}

// !!!!!!!!!!! ['__class__', '__delattr__', '__dict__', '__doc__', '__eq__', '__format__', '__getattribute__', '__hash__', '__init__', '__module__', '__ne__', '__new__', '__reduce__', '__reduce_ex__', '__repr__', '__setattr__', '__sizeof__', '__str__', '__subclasshook__', '__weakref__', '_parsed_fields', 'constants', 'fields', 'full_name', 'has_header', 'header_present', 'names', 'package', 'parsed_fields', 'short_name', 'text', 'types']




template <class ContainerAllocator>
struct IsFixedSize< ::wam_msgs::RTOrtnVel_<ContainerAllocator> >
  : TrueType
  { };

template <class ContainerAllocator>
struct IsFixedSize< ::wam_msgs::RTOrtnVel_<ContainerAllocator> const>
  : TrueType
  { };

template <class ContainerAllocator>
struct IsMessage< ::wam_msgs::RTOrtnVel_<ContainerAllocator> >
  : TrueType
  { };

template <class ContainerAllocator>
struct IsMessage< ::wam_msgs::RTOrtnVel_<ContainerAllocator> const>
  : TrueType
  { };

template <class ContainerAllocator>
struct HasHeader< ::wam_msgs::RTOrtnVel_<ContainerAllocator> >
  : FalseType
  { };

template <class ContainerAllocator>
struct HasHeader< ::wam_msgs::RTOrtnVel_<ContainerAllocator> const>
  : FalseType
  { };


template<class ContainerAllocator>
struct MD5Sum< ::wam_msgs::RTOrtnVel_<ContainerAllocator> >
{
  static const char* value()
  {
    return "2326f85574083a0a1fc4fddeff59781c";
  }

  static const char* value(const ::wam_msgs::RTOrtnVel_<ContainerAllocator>&) { return value(); }
  static const uint64_t static_value1 = 0x2326f85574083a0aULL;
  static const uint64_t static_value2 = 0x1fc4fddeff59781cULL;
};

template<class ContainerAllocator>
struct DataType< ::wam_msgs::RTOrtnVel_<ContainerAllocator> >
{
  static const char* value()
  {
    return "wam_msgs/RTOrtnVel";
  }

  static const char* value(const ::wam_msgs::RTOrtnVel_<ContainerAllocator>&) { return value(); }
};

template<class ContainerAllocator>
struct Definition< ::wam_msgs::RTOrtnVel_<ContainerAllocator> >
{
  static const char* value()
  {
    return "float32[3] angular\n\
float32	   magnitude\n\
\n\
";
  }

  static const char* value(const ::wam_msgs::RTOrtnVel_<ContainerAllocator>&) { return value(); }
};

} // namespace message_traits
} // namespace ros

namespace ros
{
namespace serialization
{

  template<class ContainerAllocator> struct Serializer< ::wam_msgs::RTOrtnVel_<ContainerAllocator> >
  {
    template<typename Stream, typename T> inline static void allInOne(Stream& stream, T m)
    {
      stream.next(m.angular);
      stream.next(m.magnitude);
    }

    ROS_DECLARE_ALLINONE_SERIALIZER;
  }; // struct RTOrtnVel_

} // namespace serialization
} // namespace ros

namespace ros
{
namespace message_operations
{

template<class ContainerAllocator>
struct Printer< ::wam_msgs::RTOrtnVel_<ContainerAllocator> >
{
  template<typename Stream> static void stream(Stream& s, const std::string& indent, const ::wam_msgs::RTOrtnVel_<ContainerAllocator>& v)
  {
    s << indent << "angular[]" << std::endl;
    for (size_t i = 0; i < v.angular.size(); ++i)
    {
      s << indent << "  angular[" << i << "]: ";
      Printer<float>::stream(s, indent + "  ", v.angular[i]);
    }
    s << indent << "magnitude: ";
    Printer<float>::stream(s, indent + "  ", v.magnitude);
  }
};

} // namespace message_operations
} // namespace ros

#endif // WAM_MSGS_MESSAGE_RTORTNVEL_H
