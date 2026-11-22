#ifndef XGALOIS_TEST_ORACLE_TINY_JSON_HPP_
#define XGALOIS_TEST_ORACLE_TINY_JSON_HPP_

#include <charconv>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <map>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <utility>
#include <variant>
#include <vector>

namespace xg {
namespace oracle {

class JsonValue {
 public:
  using Array = std::vector<JsonValue>;
  using Object = std::map<std::string, JsonValue>;
  using Storage = std::variant<std::nullptr_t, bool, int64_t, double,
                               std::string, Array, Object>;

  JsonValue() noexcept : data_(nullptr) {}
  JsonValue(std::nullptr_t) noexcept : data_(nullptr) {}
  JsonValue(bool value) noexcept : data_(value) {}
  JsonValue(int value) noexcept : data_(static_cast<int64_t>(value)) {}
  JsonValue(int64_t value) noexcept : data_(value) {}
  JsonValue(double value) noexcept : data_(value) {}
  JsonValue(const char* value) : data_(std::string(value)) {}
  JsonValue(std::string_view value) : data_(std::string(value)) {}
  JsonValue(std::string value) noexcept : data_(std::move(value)) {}
  JsonValue(Array value) noexcept : data_(std::move(value)) {}
  JsonValue(Object value) noexcept : data_(std::move(value)) {}

  bool IsNull() const noexcept {
    return std::holds_alternative<std::nullptr_t>(data_);
  }
  bool IsBool() const noexcept { return std::holds_alternative<bool>(data_); }
  bool IsInt() const noexcept { return std::holds_alternative<int64_t>(data_); }
  bool IsNumber() const noexcept {
    return std::holds_alternative<int64_t>(data_) ||
           std::holds_alternative<double>(data_);
  }
  bool IsString() const noexcept {
    return std::holds_alternative<std::string>(data_);
  }
  bool IsArray() const noexcept { return std::holds_alternative<Array>(data_); }
  bool IsObject() const noexcept {
    return std::holds_alternative<Object>(data_);
  }

  bool AsBool() const {
    if (!IsBool()) {
      throw std::runtime_error("JsonValue::AsBool: value is not a boolean");
    }
    return std::get<bool>(data_);
  }

  int64_t AsInt64() const {
    if (!IsInt()) {
      throw std::runtime_error("JsonValue::AsInt64: value is not an integer");
    }
    return std::get<int64_t>(data_);
  }

  uint64_t AsUint64() const {
    const int64_t v = AsInt64();
    if (v < 0) {
      throw std::runtime_error("JsonValue::AsUint64: integer is negative");
    }
    return static_cast<uint64_t>(v);
  }

  double AsDouble() const {
    if (IsInt()) {
      return static_cast<double>(std::get<int64_t>(data_));
    }
    if (std::holds_alternative<double>(data_)) {
      return std::get<double>(data_);
    }
    throw std::runtime_error("JsonValue::AsDouble: value is not a number");
  }

  const std::string& AsString() const {
    if (!IsString()) {
      throw std::runtime_error("JsonValue::AsString: value is not a string");
    }
    return std::get<std::string>(data_);
  }

  const Array& AsArray() const {
    if (!IsArray()) {
      throw std::runtime_error("JsonValue::AsArray: value is not an array");
    }
    return std::get<Array>(data_);
  }

  const Object& AsObject() const {
    if (!IsObject()) {
      throw std::runtime_error("JsonValue::AsObject: value is not an object");
    }
    return std::get<Object>(data_);
  }

  const JsonValue& At(const std::string& key) const {
    const Object& obj = AsObject();
    const auto it = obj.find(key);
    if (it == obj.end()) {
      throw std::runtime_error("JsonValue::At: missing key '" + key + "'");
    }
    return it->second;
  }

  const JsonValue& At(std::size_t index) const {
    const Array& arr = AsArray();
    if (index >= arr.size()) {
      throw std::runtime_error("JsonValue::At: array index out of bounds");
    }
    return arr[index];
  }

  bool Contains(const std::string& key) const {
    if (!IsObject()) {
      return false;
    }
    const Object& obj = std::get<Object>(data_);
    return obj.find(key) != obj.end();
  }

  std::size_t Size() const {
    if (IsArray()) {
      return std::get<Array>(data_).size();
    }
    if (IsObject()) {
      return std::get<Object>(data_).size();
    }
    throw std::runtime_error(
        "JsonValue::Size: value is neither an array nor an object");
  }

  bool operator==(const JsonValue& other) const = default;

 private:
  Storage data_;
};

namespace detail {

inline constexpr int kMaxParseDepth = 256;

inline void AppendUtf8(std::string& out, uint32_t cp) {
  if (cp <= 0x7FU) {
    out.push_back(static_cast<char>(cp));
  } else if (cp <= 0x7FFU) {
    out.push_back(static_cast<char>(0xC0U | ((cp >> 6) & 0x1FU)));
    out.push_back(static_cast<char>(0x80U | (cp & 0x3FU)));
  } else if (cp <= 0xFFFFU) {
    out.push_back(static_cast<char>(0xE0U | ((cp >> 12) & 0x0FU)));
    out.push_back(static_cast<char>(0x80U | ((cp >> 6) & 0x3FU)));
    out.push_back(static_cast<char>(0x80U | (cp & 0x3FU)));
  } else {
    out.push_back(static_cast<char>(0xF0U | ((cp >> 18) & 0x07U)));
    out.push_back(static_cast<char>(0x80U | ((cp >> 12) & 0x3FU)));
    out.push_back(static_cast<char>(0x80U | ((cp >> 6) & 0x3FU)));
    out.push_back(static_cast<char>(0x80U | (cp & 0x3FU)));
  }
}

class JsonParser {
 public:
  explicit JsonParser(std::string_view input) noexcept
      : input_(input), pos_(0) {}

  JsonValue ParseRoot() {
    SkipWhitespace();
    JsonValue val = ParseValue(0);
    SkipWhitespace();
    if (pos_ != input_.size()) {
      Fail("unexpected trailing characters");
    }
    return val;
  }

 private:
  [[noreturn]] void Fail(const std::string& msg) const {
    throw std::runtime_error("JSON parse error at offset " +
                             std::to_string(pos_) + ": " + msg);
  }

  bool AtEnd() const noexcept { return pos_ >= input_.size(); }

  char Peek() const {
    if (AtEnd()) {
      Fail("unexpected end of input");
    }
    return input_[pos_];
  }

  char Advance() {
    const char c = Peek();
    ++pos_;
    return c;
  }

  void SkipWhitespace() noexcept {
    while (pos_ < input_.size()) {
      const char c = input_[pos_];
      if (c == ' ' || c == '\t' || c == '\r' || c == '\n') {
        ++pos_;
      } else {
        break;
      }
    }
  }

  JsonValue ParseValue(int depth) {
    if (depth > kMaxParseDepth) {
      Fail("maximum nesting depth exceeded");
    }
    SkipWhitespace();
    const char c = Peek();
    if (c == 'n') {
      return ParseLiteral("null", JsonValue(nullptr));
    }
    if (c == 't') {
      return ParseLiteral("true", JsonValue(true));
    }
    if (c == 'f') {
      return ParseLiteral("false", JsonValue(false));
    }
    if (c == '"') {
      return JsonValue(ParseString());
    }
    if (c == '[') {
      return ParseArray(depth + 1);
    }
    if (c == '{') {
      return ParseObject(depth + 1);
    }
    if (c == '-' || (c >= '0' && c <= '9')) {
      return ParseNumber();
    }
    Fail(std::string("unexpected character '") + c + "'");
  }

  JsonValue ParseLiteral(std::string_view lit, JsonValue val) {
    if (input_.substr(pos_, lit.size()) != lit) {
      Fail("expected literal '" + std::string(lit) + "'");
    }
    pos_ += lit.size();
    return val;
  }

  uint32_t ParseHex4() {
    if (pos_ + 4 > input_.size()) {
      Fail("truncated \\uXXXX escape sequence");
    }
    uint32_t value = 0;
    for (int i = 0; i < 4; ++i) {
      const char c = input_[pos_++];
      uint32_t digit = 0;
      if (c >= '0' && c <= '9') {
        digit = static_cast<uint32_t>(c - '0');
      } else if (c >= 'a' && c <= 'f') {
        digit = static_cast<uint32_t>(10 + (c - 'a'));
      } else if (c >= 'A' && c <= 'F') {
        digit = static_cast<uint32_t>(10 + (c - 'A'));
      } else {
        Fail("invalid hex digit in \\uXXXX escape");
      }
      value = (value << 4) | digit;
    }
    return value;
  }

  std::string ParseString() {
    if (Advance() != '"') {
      Fail("expected '\"'");
    }
    std::string out;
    while (!AtEnd()) {
      const char c = Advance();
      if (c == '"') {
        return out;
      }
      if (static_cast<unsigned char>(c) < 0x20U) {
        Fail("unescaped control character in string");
      }
      if (c != '\\') {
        out.push_back(c);
        continue;
      }
      const char esc = Advance();
      switch (esc) {
        case '"':
          out.push_back('"');
          break;
        case '\\':
          out.push_back('\\');
          break;
        case '/':
          out.push_back('/');
          break;
        case 'b':
          out.push_back('\b');
          break;
        case 'f':
          out.push_back('\f');
          break;
        case 'n':
          out.push_back('\n');
          break;
        case 'r':
          out.push_back('\r');
          break;
        case 't':
          out.push_back('\t');
          break;
        case 'u': {
          uint32_t cp = ParseHex4();
          if (cp >= 0xD800U && cp <= 0xDBFFU) {
            if (pos_ + 2 > input_.size() || input_[pos_] != '\\' ||
                input_[pos_ + 1] != 'u') {
              Fail("expected low surrogate after high surrogate");
            }
            pos_ += 2;
            const uint32_t low = ParseHex4();
            if (low < 0xDC00U || low > 0xDFFFU) {
              Fail("invalid low surrogate code unit");
            }
            cp = 0x10000U + (((cp - 0xD800U) << 10) | (low - 0xDC00U));
          } else if (cp >= 0xDC00U && cp <= 0xDFFFU) {
            Fail("unexpected unpaired low surrogate");
          }
          AppendUtf8(out, cp);
          break;
        }
        default:
          Fail(std::string("invalid escape sequence '\\") + esc + "'");
      }
    }
    Fail("unterminated string literal");
  }

  JsonValue ParseNumber() {
    const std::size_t start = pos_;
    if (input_[pos_] == '-') {
      ++pos_;
    }
    if (AtEnd()) {
      Fail("incomplete number");
    }
    if (input_[pos_] == '0') {
      ++pos_;
      if (!AtEnd() && input_[pos_] >= '0' && input_[pos_] <= '9') {
        Fail("leading zeros are not allowed in JSON numbers");
      }
    } else if (input_[pos_] >= '1' && input_[pos_] <= '9') {
      while (!AtEnd() && input_[pos_] >= '0' && input_[pos_] <= '9') {
        ++pos_;
      }
    } else {
      Fail("invalid number syntax");
    }

    bool is_float = false;
    if (!AtEnd() && input_[pos_] == '.') {
      is_float = true;
      ++pos_;
      if (AtEnd() || input_[pos_] < '0' || input_[pos_] > '9') {
        Fail("expected digit after decimal point");
      }
      while (!AtEnd() && input_[pos_] >= '0' && input_[pos_] <= '9') {
        ++pos_;
      }
    }
    if (!AtEnd() && (input_[pos_] == 'e' || input_[pos_] == 'E')) {
      is_float = true;
      ++pos_;
      if (!AtEnd() && (input_[pos_] == '+' || input_[pos_] == '-')) {
        ++pos_;
      }
      if (AtEnd() || input_[pos_] < '0' || input_[pos_] > '9') {
        Fail("expected digit in exponent");
      }
      while (!AtEnd() && input_[pos_] >= '0' && input_[pos_] <= '9') {
        ++pos_;
      }
    }

    const std::string_view token = input_.substr(start, pos_ - start);
    if (!is_float) {
      int64_t int_val = 0;
      const auto [ptr, ec] =
          std::from_chars(token.data(), token.data() + token.size(), int_val);
      if (ec == std::errc{} && ptr == token.data() + token.size()) {
        return JsonValue(int_val);
      }
    }

    const std::string copy(token);
    char* end_ptr = nullptr;
    const double dbl_val = std::strtod(copy.c_str(), &end_ptr);
    if (end_ptr != copy.c_str() + copy.size() || !std::isfinite(dbl_val)) {
      Fail("number out of range: " + copy);
    }
    return JsonValue(dbl_val);
  }

  JsonValue ParseArray(int depth) {
    Advance();  // '['
    JsonValue::Array arr;
    SkipWhitespace();
    if (!AtEnd() && Peek() == ']') {
      ++pos_;
      return JsonValue(std::move(arr));
    }
    while (true) {
      arr.push_back(ParseValue(depth));
      SkipWhitespace();
      const char c = Advance();
      if (c == ']') {
        return JsonValue(std::move(arr));
      }
      if (c != ',') {
        Fail("expected ',' or ']' in array");
      }
    }
  }

  JsonValue ParseObject(int depth) {
    Advance();  // '{'
    JsonValue::Object obj;
    SkipWhitespace();
    if (!AtEnd() && Peek() == '}') {
      ++pos_;
      return JsonValue(std::move(obj));
    }
    while (true) {
      SkipWhitespace();
      if (Peek() != '"') {
        Fail("expected string key in object");
      }
      std::string key = ParseString();
      SkipWhitespace();
      if (Advance() != ':') {
        Fail("expected ':' after object key");
      }
      JsonValue val = ParseValue(depth);
      obj.insert_or_assign(std::move(key), std::move(val));
      SkipWhitespace();
      const char c = Advance();
      if (c == '}') {
        return JsonValue(std::move(obj));
      }
      if (c != ',') {
        Fail("expected ',' or '}' in object");
      }
    }
  }

  std::string_view input_;
  std::size_t pos_;
};

inline void EscapeJsonString(std::string& out, const std::string& s) {
  static constexpr char kHex[] = "0123456789abcdef";
  out.push_back('"');
  for (unsigned char c : s) {
    switch (c) {
      case '"':
        out += "\\\"";
        break;
      case '\\':
        out += "\\\\";
        break;
      case '\b':
        out += "\\b";
        break;
      case '\f':
        out += "\\f";
        break;
      case '\n':
        out += "\\n";
        break;
      case '\r':
        out += "\\r";
        break;
      case '\t':
        out += "\\t";
        break;
      default:
        if (c < 0x20U) {
          out += "\\u00";
          out.push_back(kHex[(c >> 4) & 0x0FU]);
          out.push_back(kHex[c & 0x0FU]);
        } else {
          out.push_back(static_cast<char>(c));
        }
        break;
    }
  }
  out.push_back('"');
}

inline void DumpJsonImpl(const JsonValue& v, int indent, int level,
                         std::string& out) {
  if (v.IsNull()) {
    out += "null";
  } else if (v.IsBool()) {
    out += v.AsBool() ? "true" : "false";
  } else if (v.IsInt()) {
    out += std::to_string(v.AsInt64());
  } else if (v.IsNumber()) {
    const double d = v.AsDouble();
    if (!std::isfinite(d)) {
      throw std::runtime_error(
          "DumpJson: non-finite floating-point value cannot be serialized");
    }
    std::ostringstream oss;
    oss << std::setprecision(17) << d;
    std::string s = oss.str();
    if (s.find('.') == std::string::npos && s.find('e') == std::string::npos &&
        s.find('E') == std::string::npos) {
      s += ".0";
    }
    out += s;
  } else if (v.IsString()) {
    EscapeJsonString(out, v.AsString());
  } else if (v.IsArray()) {
    const auto& arr = v.AsArray();
    if (arr.empty()) {
      out += "[]";
      return;
    }
    out.push_back('[');
    for (std::size_t i = 0; i < arr.size(); ++i) {
      if (i > 0) {
        out.push_back(',');
      }
      if (indent > 0) {
        out.push_back('\n');
        out.append(static_cast<std::size_t>((level + 1) * indent), ' ');
      }
      DumpJsonImpl(arr[i], indent, level + 1, out);
    }
    if (indent > 0) {
      out.push_back('\n');
      out.append(static_cast<std::size_t>(level * indent), ' ');
    }
    out.push_back(']');
  } else if (v.IsObject()) {
    const auto& obj = v.AsObject();
    if (obj.empty()) {
      out += "{}";
      return;
    }
    out.push_back('{');
    bool first = true;
    for (const auto& [k, val] : obj) {
      if (!first) {
        out.push_back(',');
      }
      first = false;
      if (indent > 0) {
        out.push_back('\n');
        out.append(static_cast<std::size_t>((level + 1) * indent), ' ');
      }
      EscapeJsonString(out, k);
      out += (indent > 0) ? ": " : ":";
      DumpJsonImpl(val, indent, level + 1, out);
    }
    if (indent > 0) {
      out.push_back('\n');
      out.append(static_cast<std::size_t>(level * indent), ' ');
    }
    out.push_back('}');
  }
}

}  // namespace detail

inline JsonValue ParseJson(std::string_view text) {
  detail::JsonParser parser(text);
  return parser.ParseRoot();
}

inline JsonValue ParseJsonFile(const std::string& path) {
  std::ifstream file(path, std::ios::in | std::ios::binary);
  if (!file.is_open()) {
    throw std::runtime_error("ParseJsonFile: failed to open '" + path + "'");
  }
  std::ostringstream oss;
  oss << file.rdbuf();
  if (file.bad()) {
    throw std::runtime_error("ParseJsonFile: failed to read '" + path + "'");
  }
  return ParseJson(oss.str());
}

inline std::string DumpJson(const JsonValue& v, int indent = 0) {
  std::string out;
  detail::DumpJsonImpl(v, indent, 0, out);
  return out;
}

}  // namespace oracle
}  // namespace xg

#endif  // XGALOIS_TEST_ORACLE_TINY_JSON_HPP_
