# -*- coding: utf-8 -*-
# Generated by the protocol buffer compiler.  DO NOT EDIT!
# source: narupa/protocol/command/command_service.proto
"""Generated protocol buffer code."""
from google.protobuf import descriptor as _descriptor
from google.protobuf import descriptor_pool as _descriptor_pool
from google.protobuf import symbol_database as _symbol_database
from google.protobuf.internal import builder as _builder
# @@protoc_insertion_point(imports)

_sym_db = _symbol_database.Default()


from google.protobuf import struct_pb2 as google_dot_protobuf_dot_struct__pb2


DESCRIPTOR = _descriptor_pool.Default().AddSerializedFile(b'\n-narupa/protocol/command/command_service.proto\x12\x17narupa.protocol.command\x1a\x1cgoogle/protobuf/struct.proto\"\x14\n\x12GetCommandsRequest\"M\n\x10GetCommandsReply\x12\x39\n\x08\x63ommands\x18\x01 \x03(\x0b\x32\'.narupa.protocol.command.CommandMessage\"7\n\x0c\x43ommandReply\x12\'\n\x06result\x18\x01 \x01(\x0b\x32\x17.google.protobuf.Struct\"J\n\x0e\x43ommandMessage\x12\x0c\n\x04name\x18\x01 \x01(\t\x12*\n\targuments\x18\x02 \x01(\x0b\x32\x17.google.protobuf.Struct2\xd2\x01\n\x07\x43ommand\x12g\n\x0bGetCommands\x12+.narupa.protocol.command.GetCommandsRequest\x1a).narupa.protocol.command.GetCommandsReply\"\x00\x12^\n\nRunCommand\x12\'.narupa.protocol.command.CommandMessage\x1a%.narupa.protocol.command.CommandReply\"\x00\x62\x06proto3')

_globals = globals()
_builder.BuildMessageAndEnumDescriptors(DESCRIPTOR, _globals)
_builder.BuildTopDescriptorsAndMessages(DESCRIPTOR, 'narupa.protocol.command.command_service_pb2', _globals)
if _descriptor._USE_C_DESCRIPTORS == False:
  DESCRIPTOR._options = None
  _globals['_GETCOMMANDSREQUEST']._serialized_start=104
  _globals['_GETCOMMANDSREQUEST']._serialized_end=124
  _globals['_GETCOMMANDSREPLY']._serialized_start=126
  _globals['_GETCOMMANDSREPLY']._serialized_end=203
  _globals['_COMMANDREPLY']._serialized_start=205
  _globals['_COMMANDREPLY']._serialized_end=260
  _globals['_COMMANDMESSAGE']._serialized_start=262
  _globals['_COMMANDMESSAGE']._serialized_end=336
  _globals['_COMMAND']._serialized_start=339
  _globals['_COMMAND']._serialized_end=549
# @@protoc_insertion_point(module_scope)
