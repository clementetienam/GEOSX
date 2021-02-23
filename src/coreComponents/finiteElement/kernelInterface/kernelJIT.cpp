#include "KernelBase.hpp"
#include "kernelJITCompileCommands.hpp"

namespace geosx
{
namespace finiteElement
{

jitti::CompilationInfo getCompilationInfo()
{
  jitti::CompilationInfo info;

  info.compileCommand = kernelJIT_COMPILE_COMMAND;
  info.linker = kernelJIT_LINKER;
  info.linkArgs = kernelJIT_LINK_ARGS;

  std::string const currentFile = __FILE__;
  info.header = currentFile.substr( 0, currentFile.size() - ( sizeof( "kernelJIT.cpp" ) - 1 ) )
                                 + "KernelBase.hpp";

  info.function = "geosx::finiteElement::fooBar";

  return info;
}

} // namespace finiteElement
} // namespace geosx
