from mdvtools.dbutils.mdv_server_app import app

"""
As of this writing, this is a way of running a parallel instance of the server inside a debugger (debug this file from inside container).

Attempting to view a project in Safari, I get something like this:

127.0.0.1 - - [10/Feb/2025 16:16:48] code 400, message Bad request version ('³`ØP\x9e=')
127.0.0.1 - - [10/Feb/2025 16:16:48] "\x16\x03\x01\x02\x00\x01\x00\x01ü\x03\x03\x94Ui£ê\x97\x0c°,T\x1cQ>Ê\x9fiÜ¾kMÇF\x1b_úF)¼\x194}å \x11cÓ]Ça\x1f³`ØP\x9e=" 400 -
127.0.0.1 - - [10/Feb/2025 16:16:48] code 400, message Bad HTTP/0.9 request type ('\x16\x03\x01\x02\x00\x01\x00\x01ü\x03\x03@½á\x97¹`\x03HMNÍôã\x9cKW#üï\x9e1-\x87.GZ~\x81Ò[\x8b\x81')
127.0.0.1 - - [10/Feb/2025 16:16:48] "\x16\x03\x01\x02\x00\x01\x00\x01ü\x03\x03@½á\x97¹`\x03HMNÍôã\x9cKW#üï\x9e1-\x87.GZ~\x81Ò[\x8b\x81 ñ&M"Ri\x13\x875Z£\x1b¯\x18Û.#s[Ë¼ýÓzO>jïV\x03ÉÄ\x00*ÊÊ\x13\x01\x13\x02\x13\x03À,À+Ì©À0À/Ì¨À" 400 -
127.0.0.1 - - [10/Feb/2025 16:16:48] code 400, message Bad request version ('¥`®_sf1Îä\x00*ªª\x13\x01\x13\x02\x13\x03À,À+Ì©À0À/Ì¨À')
127.0.0.1 - - [10/Feb/2025 16:16:48] "\x16\x03\x01\x02\x00\x01\x00\x01ü\x03\x03\x9fñG8<\x0fû]@oú\x91\x83¡@\x1bø8KX\x8e8!\x9a\x04åÉ´\x09|\x19_ \x02G\x87Ì\x18\x8fK%$@«ç\x99Ôd\x8d{\x98}ÛI¸ ¥`®_sf1Îä\x00*ªª\x13\x01\x13\x02\x13\x03À,À+Ì©À0À/Ì¨À" 400 -
127.0.0.1 - - [10/Feb/2025 16:16:48] code 400, message Bad request version ('\x17È¸c®ÍêTw\x9b\x00*ªª\x13\x01\x13\x02\x13\x03À,À+Ì©À0À/Ì¨À')
127.0.0.1 - - [10/Feb/2025 16:16:48] "\x16\x03\x01\x02\x00\x01\x00\x01ü\x03\x03÷%Nc\x91\x82IØ§W\x07ç\x88`|\x0dÃo-\x05O\x12=wo1\x8eT\x86\x91ä\x9a wY\x95³ã"³@±¡Éï+&8~È_W\x96\x1a\x1e\x17È¸c®ÍêTw\x9b\x00*ªª\x13\x01\x13\x02\x13\x03À,À+Ì©À0À/Ì¨À" 400 -
127.0.0.1 - - [10/Feb/2025 16:16:48] code 400, message Bad request version ('\x08t\x04òã`\x8ee\x07j/r]\x13â¯¯z6R½\x88õÿ\x01ªÎ5\x00*\x1a\x1a\x13\x01\x13\x02\x13\x03À,À+Ì©À0À/Ì¨À')
127.0.0.1 - - [10/Feb/2025 16:16:48] "\x16\x03\x01\x02\x00\x01\x00\x01ü\x03\x03\x9a§CUÝ\x12Fhë\x05ú(\x13Ï5oìönÚ_\x07:;\x1f)\x1f\x97ÒÔ\x1b´ lþ\x0b \x08t\x04òã`\x8ee\x07j/r]\x13â¯¯z6R½\x88õÿ\x01ªÎ5\x00*\x1a\x1a\x13\x01\x13\x02\x13\x03À,À+Ì©À0À/Ì¨À" 400 -
127.0.0.1 - - [10/Feb/2025 16:16:48] code 400, message Bad HTTP/0.9 request type ('\x16\x03\x01\x02\x00\x01\x00\x01ü\x03\x03\x9bøkd\x89ÅêõÞï¯\x95u%\x1b*wu¶ð\x1a$\x8d9µôk6\x97ë\x9fX')
127.0.0.1 - - [10/Feb/2025 16:16:48] "\x16\x03\x01\x02\x00\x01\x00\x01ü\x03\x03\x9bøkd\x89ÅêõÞï¯\x95u%\x1b*wu¶ð\x1a$\x8d9µôk6\x97ë\x9fX È\x0eT\x9c¡´é_Ò2/Ii\x19\x86\x8c|»#\x92ÒÛ®7²\x8a¶EvÀ§}\x00*ZZ\x13\x01\x13\x02\x13\x03À,À+Ì©À0À/Ì¨À" 400 -
127.0.0.1 - - [10/Feb/2025 16:16:48] code 400, message Bad request version ('º\x99fO^Y\x95\x95Ñ»Î\x83ú¢\x9f¬ìÎá¹\x8bõ\x8a;LfÐ\x9f$&ôî\x00*::\x13\x01\x13\x02\x13\x03À,À+Ì©À0À/Ì¨À')
127.0.0.1 - - [10/Feb/2025 16:16:48] "\x16\x03\x01\x02\x00\x01\x00\x01ü\x03\x03mI\x9e\x01\x01²³\x88ñ½º]ù\x83\x91  SÛ\x8diê¾ï9^Ù[\x87¹<y º\x99fO^Y\x95\x95Ñ»Î\x83ú¢\x9f¬ìÎá¹\x8bõ\x8a;LfÐ\x9f$&ôî\x00*::\x13\x01\x13\x02\x13\x03À,À+Ì©À0À/Ì¨À" 400 -
127.0.0.1 - - [10/Feb/2025 16:16:49] code 400, message Bad request version ('\x0fp\x00*::\x13\x01\x13\x02\x13\x03À,À+Ì©À0À/Ì¨À')
127.0.0.1 - - [10/Feb/2025 16:16:49] "\x16\x03\x01\x02\x00\x01\x00\x01ü\x03\x03w¿<ÿËQÝêa\x93,RµýævtÄä\x97äº\x0bp_"nÒJ\x87k^ \x1fçíÛ\x07ì¿á¹\x91\x89o³\x7fRE\x85Ã(ë7\x0bR»R0\x90¶å\x0c\x0fp\x00*::\x13\x01\x13\x02\x13\x03À,À+Ì©À0À/Ì¨À" 400 -

Couldn't see much by stopping on a breakpoint and examining the request object.  I'm not sure what's going on here.
But it's good to have a way to run the server in a debugger.
"""


if __name__ == '__main__':
    app.run(debug=True, port=5056)