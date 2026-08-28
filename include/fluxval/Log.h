// Log.h -- print to the console and record the same text for the .txt summary,
// so the file in the output directory is exactly what you saw on screen.

#pragma once

#include <string>
#include <vector>

namespace fluxval {

std::vector<std::string>& logLines();

/// printf to stdout and append the formatted line to the record.
void logf(const char* fmt, ...);

/// Append a line that is already formatted.
void logLine(const std::string& s);

/// Append without printing (for headers that only belong in the file).
void logQuiet(const std::string& s);

/// Write everything recorded so far to `path`.
bool writeLog(const std::string& path);

}  // namespace fluxval
