"""
PyPair: DNA Sequence Encoding/Decoding Module
==============================================

Provides 2-bit encoding for DNA sequences and run-length encoding
for compression and feature extraction.

Usage:
    from models.encoding import Encode_decode

    encoder = Encode_decode()
    encoded = encoder.encode("ACGTACGT")  # "0001101100010110"
    decoded = encoder.decode(encoded)      # "ACGTACGT"
"""

from collections import OrderedDict


class Encode_decode:
    """
    Encoder/Decoder class for 2-bit encoding of DNA sequences.

    2-bit encoding reduces storage requirements:
    - DNA: 4 characters (A, C, G, T)
    - 2 bits can represent 4 states (00, 01, 10, 11)
    - 28bp k-mer: 28 bytes → 56 bits (7 bytes)
    """

    # Encoding dictionary: DNA → 2-bit binary
    ENCODE_DICT = {
        '$': '',     # Sentinel character (removed)
        ',': '',     # Comma (removed)
        'A': '00',   # Adenine
        'C': '01',   # Cytosine
        'G': '10',   # Guanine
        'T': '11'    # Thymine
    }

    # Decoding dictionary: 2-bit binary → DNA
    DECODE_DICT = {
        '00': 'A',
        '01': 'C',
        '10': 'G',
        '11': 'T'
    }

    def encode(self, index):
        """
        Encode DNA sequence(s) using 2-bit encoding.

        Args:
            index: String or list of strings to encode

        Returns:
            Encoded string or list of encoded strings

        Example:
            >>> encoder = Encode_decode()
            >>> encoder.encode("ACGT")
            '00011011'
            >>> encoder.encode(["ACG", "TGC"])
            ['000110', '110110']
        """
        encoded = []

        if isinstance(index, list):
            # Encode list of strings
            for substring in index:
                trans_table = substring.maketrans(self.ENCODE_DICT)
                txt = substring.translate(trans_table)
                encoded.append(txt)
        else:
            # Encode single string
            trans_table = index.maketrans(self.ENCODE_DICT)
            txt = index.translate(trans_table)
            encoded = txt

        return encoded

    def decode(self, index):
        """
        Decode 2-bit encoded sequence(s) back to DNA.

        Args:
            index: Encoded string or list of encoded strings

        Returns:
            Decoded DNA string or list of DNA strings

        Example:
            >>> encoder = Encode_decode()
            >>> encoder.decode("00011011")
            'ACGT'
            >>> encoder.decode(["000110", "110110"])
            ['ACG', 'TGC']
        """
        decoded = []
        decoded_substring = []

        if isinstance(index, list):
            # Decode list of strings
            for i in range(len(index)):
                count = 0
                while count < len(index):
                    substring = index[count] + index[count + 1]
                    decoded_substring.append(
                        self.DECODE_DICT.get(str(substring), str(substring))
                    )
                    count += 2
                decoded_string = ''.join(decoded_substring)
                decoded.append(decoded_string)
                decoded_substring = []  # Reset for next iteration
        else:
            # Decode single string
            count = 0
            while count < len(index):
                substring = index[count] + index[count + 1]
                decoded_substring.append(
                    self.DECODE_DICT.get(str(substring), str(substring))
                )
                count += 2
            decoded_string = ''.join(decoded_substring)
            decoded = decoded_string

        return decoded

    @staticmethod
    def run_length_encode(input_string):
        """
        Run-length encoding of DNA sequence for compression.

        Encodes consecutive identical characters as character + count.

        Args:
            input_string: DNA sequence string

        Returns:
            tuple: (output, characters, values)
                - output: Full RLE string (e.g., "A3C2G1")
                - characters: Just the characters (e.g., "ACG")
                - values: Just the counts (e.g., "321")

        Example:
            >>> Encode_decode.run_length_encode("AAACGGG")
            ('A3C1G3', 'ACG', '313')
        """
        count_dict = OrderedDict.fromkeys(input_string, 0)

        for ch in input_string:
            count_dict[ch] += 1

        output = ''
        characters = ''
        values = ''

        for key, value in count_dict.items():
            output += key + str(value)
            characters += key
            values += str(value)

        return output, characters, values


# Backwards compatibility - keep old function name
def runLengthEncoding(input_string):
    """
    Legacy function name for run-length encoding.
    Use Encode_decode.run_length_encode() instead.
    """
    output, characters, values = Encode_decode.run_length_encode(input_string)
    return output
