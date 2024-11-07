#ifndef ZFILEREADER_H
#define ZFILEREADER_H

// taken from web:
// https://stackoverflow.com/questions/7273326/getting-the-nth-line-of-a-text-file-in-c?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
// https://stackoverflow.com/questions/236129/the-most-elegant-way-to-iterate-the-words-of-a-string?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa
// https://stackoverflow.com/questions/8362094/replace-multiple-spaces-with-one-space-in-a-string?utm_medium=organic&utm_source=google_rich_qa&utm_campaign=google_rich_qa

#include <string>
#include <sstream>
#include <vector>
#include <iterator>

template<typename Out>
void split(const std::string &s, char delim, Out result) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        *(result++) = item;
    }
}

std::vector<std::string> split(const std::string &s, char delim);

// bool BothAreSpaces(char lhs, char rhs);

#include <fstream>

class TextFile {
    std::ifstream file_stream;
    std::vector<std::streampos> linebegins;
    TextFile& operator=(TextFile& b) = delete;
    const int some_reasonable_max_line_length = 1024;
public:
    TextFile(std::string filename, const int nMaxLines = -1)
    :file_stream(filename)
    {
        //this chunk stolen from Armen's,
        std::string s;
        //for performance
        s.reserve(some_reasonable_max_line_length);
        int l = 0;
        while(file_stream) {
            if(nMaxLines >= 0 && l >= nMaxLines)
              break;
            linebegins.push_back(file_stream.tellg());
            std::getline(file_stream, s);
            l++;
        }
    }
    TextFile(TextFile&& b) : file_stream(std::move(b.file_stream)), linebegins(std::move(b.linebegins)){
    }
    void operator=(TextFile&& b){
        file_stream = std::move(b.file_stream);
        linebegins = std::move(b.linebegins);
    }
    std::string ReadNthLine(unsigned int N) {
        if (N >= linebegins.size()-1)
            throw std::runtime_error("File doesn't have that many lines!");
        std::string s;
        // clear EOF and error flags
        file_stream.clear();
        file_stream.seekg(linebegins[N]);
        std::getline(file_stream, s);
        return s;
    }
};

#endif // ZFILEREADER_H
