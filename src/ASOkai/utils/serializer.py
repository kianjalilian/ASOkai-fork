#!/usr/bin/env python
"""
Filename: src/ASOkai/utils/serializer.py
Author: Arash Ayat
Copyright: 2025, Alexander Schliep
Version: 0.1.0
Description: This file defines the Serializable class for serializing and deserializing objects.
License: LGPL-3.0-or-later
"""
from typing import Dict, Any, Type, TypeVar, Callable
import json
import importlib
import inspect

T = TypeVar('T', bound='Serializable')


class Serializable:
    """
    Base class for objects that can be serialized to a dictionary and file.
    """
    
    # Class-level registry for external types
    # Maps type -> (type_name, serialize_func, flatten)
    _type_serializers: Dict[type, tuple[str, Callable[[Any], Dict[str, Any]], bool]] = {}
    # Maps type_name -> deserialize_func
    _type_deserializers: Dict[str, Callable[[Dict[str, Any]], Any]] = {}
    
    @classmethod
    def register_type(
        cls,
        type_class: type,
        type_name: str,
        serialize: Callable[[Any], Dict[str, Any]],
        deserialize: Callable[[Dict[str, Any]], Any],
        flatten: bool = False
    ) -> None:
        """
        Register an external type for serialization/deserialization.
        
        Args:
            type_class: The type to register (e.g., Locus, Seq).
            type_name: A unique string identifier for the type (e.g., 'Locus', 'Bio.Seq.Seq').
            serialize: A function that takes an instance and returns a dict of its data.
            deserialize: A function that takes a dict and returns an instance.
            flatten: If True, serialize components as top-level keys instead of nested dict.
                     When flattened, no __type__ marker is added (components merge into parent).
        
        Example:
            # Nested (default): {"locus": {"chr": "12", ..., "__type__": "Locus"}}
            Serializable.register_type(
                Locus, 'Locus',
                serialize=lambda l: {'chr': l.chr, 'start': l.start, 'end': l.end, 'strand': l.strand},
                deserialize=lambda d: Locus(**d)
            )
            
            # Flattened: {"chr": "12", "start": 100, "end": 200, "strand": "-"}
            Serializable.register_type(
                Locus, 'Locus',
                serialize=lambda l: {'chr': l.chr, 'start': l.start, 'end': l.end, 'strand': l.strand},
                deserialize=lambda d: Locus(**d),
                flatten=True
            )
        """
        cls._type_serializers[type_class] = (type_name, serialize, flatten)
        cls._type_deserializers[type_name] = deserialize

    def __init__(self, **kwargs: Any) -> None:
        """
        Initialize with keyword arguments.
        All kwargs are set as instance attributes.
        """
        for key, value in kwargs.items():
            setattr(self, key, value)

    # Override in subclass to exclude attributes from serialization
    _non_serializable_attrs: set = set()

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert the object to a dictionary.
        """
        data: Dict[str, Any] = {
            '__class__': self.__class__.__name__,
            '__module__': self.__class__.__module__
        }
        for key, value in self.__dict__.items():
            if key in self._non_serializable_attrs:
                continue
            # Check if this is a flattened registered type
            flattened = False
            for type_class, (type_name, serialize_func, flatten) in self._type_serializers.items():
                if isinstance(value, type_class) and flatten:
                    # Merge components directly into parent dict
                    data.update(serialize_func(value))
                    flattened = True
                    break
            if not flattened:
                data[key] = self._serialize_value(value)
        return data

    def _serialize_value(self, value: Any) -> Any:
        if isinstance(value, (str, int, float, bool, type(None))):
            return value
        elif isinstance(value, Serializable):
            return value.to_dict()
        elif isinstance(value, dict):
            return {k: self._serialize_value(v) for k, v in value.items()}
        elif isinstance(value, (list, tuple)):
            return [self._serialize_value(item) for item in value]
        else:
            # Check registered external types
            for type_class, (type_name, serialize_func, flatten) in self._type_serializers.items():
                if isinstance(value, type_class):
                    data = serialize_func(value)
                    if not flatten:
                        data['__type__'] = type_name
                    return data
            raise TypeError(f"Cannot serialize type {type(value).__name__}")

    def to_file(self, file_path: str) -> None:
        """
        method to write the object to a file.
        """
        with open(file_path, 'w') as f:
            json.dump(self.to_dict(), f, indent=4)

    @classmethod
    def _get_init_arg_name_map(cls) -> Dict[str, str]:
        """
        Provides a mapping from serialized attribute names to __init__ parameter names.
        Override this in a subclass if the names differ.
        Example: return {"_sequence": "sequence"}
        """
        return {}

    @classmethod
    def _deserialize_value(cls, value_data: Any) -> Any:
        """
        Hook for subclasses to deserialize special value types.
        The base implementation handles registered types, Serializable objects, dicts, and lists.
        """
        if isinstance(value_data, list):
            return [cls._deserialize_value(item) for item in value_data]
        elif isinstance(value_data, dict):
            # Check for registered external type
            if '__type__' in value_data:
                type_name = value_data['__type__']
                if type_name in cls._type_deserializers:
                    # Create a copy without __type__ for the deserializer
                    data = {k: v for k, v in value_data.items() if k != '__type__'}
                    return cls._type_deserializers[type_name](data)
            # Check for Serializable object
            if '__class__' in value_data and '__module__' in value_data:
                return Serializable.from_dict(value_data)
            # Regular dict - recursively deserialize values
            return {k: cls._deserialize_value(v) for k, v in value_data.items()}
        return value_data

    @classmethod
    def from_dict(cls: Type[T], data: Dict[str, Any]) -> T:
        """
        Create an object from a dictionary using introspection.
        """
        if not (isinstance(data, dict) and '__class__' in data and '__module__' in data):
            raise ValueError("Data is not a valid serialized object dictionary.")

        # Copy to avoid mutating the input
        data = data.copy()
        module_name = data.pop('__module__')
        class_name = data.pop('__class__')
        module = importlib.import_module(module_name)
        cls_from_data = getattr(module, class_name)

        if not issubclass(cls_from_data, Serializable):
            raise TypeError(f"Class {class_name} is not a subclass of Serializable")

        if not issubclass(cls_from_data, cls):
            raise TypeError(f"Data is for class {class_name}, which is not a subclass of {cls.__name__}")

        sig = inspect.signature(cls_from_data.__init__)
        name_map = cls_from_data._get_init_arg_name_map()
        init_args: Dict[str, Any] = {}

        # First, handle mapped arguments
        for attr_name, param_name in name_map.items():
            if attr_name in data and param_name in sig.parameters:
                value = data.pop(attr_name)
                init_args[param_name] = cls_from_data._deserialize_value(value)

        # Handle remaining arguments where attr_name == param_name
        for param_name in sig.parameters:
            if param_name in data:
                value = data.pop(param_name)
                init_args[param_name] = cls_from_data._deserialize_value(value)

        # Pass any leftover data as kwargs, if __init__ accepts them
        for param in sig.parameters.values():
            if param.kind == param.VAR_KEYWORD:
                # Deserialize remaining data for kwargs
                deserialized_kwargs = {k: cls_from_data._deserialize_value(v) for k, v in data.items()}
                init_args.update(deserialized_kwargs)
                break

        return cls_from_data(**init_args)


    @classmethod
    def from_file(cls: Type[T], file_path: str) -> T:
        """
        Create an object from a JSON file.
        """
        with open(file_path, 'r') as f:
            data = json.load(f)
        return cls.from_dict(data)
