# -*- coding: utf-8 -*-
# Generated by the protocol buffer compiler.  DO NOT EDIT!
# source: narupa/protocol/instance/connect_to_trajectory.proto
"""Generated protocol buffer code."""
from google.protobuf import descriptor as _descriptor
from google.protobuf import descriptor_pool as _descriptor_pool
from google.protobuf import symbol_database as _symbol_database
from google.protobuf.internal import builder as _builder
# @@protoc_insertion_point(imports)

_sym_db = _symbol_database.Default()


from narupa.protocol import address_pb2 as narupa_dot_protocol_dot_address__pb2


DESCRIPTOR = _descriptor_pool.Default().AddSerializedFile(b'\n4narupa/protocol/instance/connect_to_trajectory.proto\x1a\x1dnarupa/protocol/address.proto\"\\\n\x1a\x43onnectToTrajectoryRequest\x12\x13\n\x0binstance_id\x18\x01 \x01(\t\x12)\n\x07\x61\x64\x64ress\x18\x02 \x01(\x0b\x32\x18.narupa.protocol.Address\"\x1d\n\x1b\x43onnectToTrajectoryResponseb\x06proto3')

_globals = globals()
_builder.BuildMessageAndEnumDescriptors(DESCRIPTOR, _globals)
_builder.BuildTopDescriptorsAndMessages(DESCRIPTOR, 'narupa.protocol.instance.connect_to_trajectory_pb2', _globals)
if _descriptor._USE_C_DESCRIPTORS == False:
  DESCRIPTOR._options = None
  _globals['_CONNECTTOTRAJECTORYREQUEST']._serialized_start=87
  _globals['_CONNECTTOTRAJECTORYREQUEST']._serialized_end=179
  _globals['_CONNECTTOTRAJECTORYRESPONSE']._serialized_start=181
  _globals['_CONNECTTOTRAJECTORYRESPONSE']._serialized_end=210
# @@protoc_insertion_point(module_scope)