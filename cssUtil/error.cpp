
/*
    pbrt source code Copyright(c) 1998-2012 Matt Pharr and Greg Humphreys.

    This file is part of pbrt.

    Redistribution and use in source and binary forms, with or without
    modification, are permitted provided that the following conditions are
    met:

    - Redistributions of source code must retain the above copyright
      notice, this list of conditions and the following disclaimer.

    - Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.

    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS
    IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
    TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A
    PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
    HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
    LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
    DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
    THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
    (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
    OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

 */



// Error Reporting Includes
#include <stdarg.h>
#include "error.h"

// Error Reporting Definitions
#define CCS_ERROR_IGNORE 0
#define CCS_ERROR_CONTINUE 1
#define CCS_ERROR_ABORT 2

int TerminalWidth() {
#if defined(CCS_IS_WINDOWS)
	HANDLE h = GetStdHandle(STD_OUTPUT_HANDLE);
	if (h == INVALID_HANDLE_VALUE || h == NULL) {
		fprintf(stderr, "GetStdHandle() call failed");
		return 80;
	}
	CONSOLE_SCREEN_BUFFER_INFO bufferInfo = { 0 };
	GetConsoleScreenBufferInfo(h, &bufferInfo);
	return bufferInfo.dwSize.X;
#else
	struct winsize w;
	if (ioctl(STDOUT_FILENO, TIOCGWINSZ, &w) < 0) {
		fprintf(stderr, "Error in ioctl() in TerminalWidth(): %d", errno);
		return 80;
	}
	return w.ws_col;
#endif // CCS_IS_WINDOWS
}

const char *findWordEnd(const char *buf) {
    while (*buf != '\0' && !isspace(*buf))
        ++buf;
    return buf;
}

// Error Reporting Functions
static void processError(const char *format, va_list args,
        const char *errorType, int disposition) {
    // Report error
    if (disposition == CCS_ERROR_IGNORE)
        return;

    // Build up an entire formatted error string and print it all at once;
    // this way, if multiple threads are printing messages at once, they
    // don't get jumbled up...
    std::string errorString;

    // PBRT_ERROR_CONTINUE, PBRT_ERROR_ABORT
    // Print formatted error message
    int width = std::max(20, TerminalWidth() - 2);
    errorString += errorType;
    errorString += ": ";
    int column = errorString.size();

#if !defined(CCS_IS_WINDOWS)
    char *errorBuf;
    if (vasprintf(&errorBuf, format, args) == -1) {
        fprintf(stderr, "vasprintf() unable to allocate memory!\n");
        abort();
    }
#else
    char errorBuf[2048];
    vsnprintf_s(errorBuf, sizeof(errorBuf), _TRUNCATE, format, args);
#endif

    const char *msgPos = errorBuf;
    while (true) {
        while (*msgPos != '\0' && isspace(*msgPos))
            ++msgPos;
        if (*msgPos == '\0')
            break;

        const char *wordEnd = findWordEnd(msgPos);
        if (column + wordEnd - msgPos > width) {
            errorString += "\n    ";
            column = 4;
        }
        while (msgPos != wordEnd) {
            errorString += *msgPos++;
            ++column;
        }
        errorString += ' ';
        ++column;
    }

    fprintf(stderr, "%s\n", errorString.c_str());

    if (disposition == CCS_ERROR_ABORT) {
#if defined(CCS_IS_WINDOWS)
        __debugbreak();
#else
        abort();
#endif
    }
#if !defined(CCS_IS_WINDOWS)
    free(errorBuf);
#endif
}


void Info(const char *format, ...) {
    va_list args;
    va_start(args, format);
    processError(format, args, "Notice", CCS_ERROR_CONTINUE);
    va_end(args);
}


void Warning(const char *format, ...) {
    va_list args;
    va_start(args, format);
    processError(format, args, "Warning", CCS_ERROR_CONTINUE);
    va_end(args);
}


void Error(const char *format, ...) {
    va_list args;
    va_start(args, format);
    processError(format, args, "Error", CCS_ERROR_CONTINUE);
    va_end(args);
}


void Severe(const char *format, ...) {
    va_list args;
    va_start(args, format);
    processError(format, args, "Fatal Error", CCS_ERROR_ABORT);
    va_end(args);
}


