#include "fluxval/Log.h"

#include <cstdarg>
#include <cstdio>
#include <fstream>

namespace fluxval {

std::vector<std::string>& logLines() {
  static std::vector<std::string> lines;
  return lines;
}

void logLine(const std::string& s) {
  logLines().push_back(s);
  std::printf("%s\n", s.c_str());
}

void logQuiet(const std::string& s) { logLines().push_back(s); }

void logf(const char* fmt, ...) {
  char b[2048];
  va_list ap;
  va_start(ap, fmt);
  std::vsnprintf(b, sizeof(b), fmt, ap);
  va_end(ap);
  std::string s(b);
  // logf() is used with trailing newlines in places; strip them so the record
  // stays one entry per line.
  while (!s.empty() && (s.back() == '\n' || s.back() == '\r')) s.pop_back();
  logLine(s);
}

bool writeLog(const std::string& path) {
  std::ofstream out(path);
  if (!out) return false;
  for (const std::string& s : logLines()) out << s << '\n';
  return true;
}

}  // namespace fluxval
